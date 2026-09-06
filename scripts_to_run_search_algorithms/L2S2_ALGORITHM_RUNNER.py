#!/usr/bin/env python3

################# L2S2_ALGORITHM_RUNNER.py #####################
#
# perturbagen_group  ---  name of signatures perturbagen's group (compulsory option), e.g. MEK, mTOR, PI3K...
# --project_path  ---  path to a directory with the whole project
# --out_folder  ---  NAME of out-folder, where will be located obtained results of this script. the location of this folder has already been determined by the script.
#
################################################################




import pandas as pd
import requests
from requests.exceptions import HTTPError
import json
import sys
import os
import argparse

url = "https://l2s2.maayanlab.cloud/graphql"  # Внимание: HTTPS!

def get_l2s2_condition_annotations(geneset: list, first=20000):
    """
    Загружает сигнатуру в L2S2 и возвращает condition annotations,
    извлечённые из поля 'term' через парсинг.
    """
    query = {
        "operationName": "EnrichmentQuery",
        "variables": {
            "genes": geneset,
            "filterTerm": "",
            "offset": 0,
            "first": first,
            "filterFda": False,
            "sortBy": "pvalue_up",
            "filterKo": False,
        },
        "query": """query EnrichmentQuery(
            $genes: [String]!,
            $filterTerm: String = "",
            $offset: Int = 0,
            $first: Int = 10,
            $filterFda: Boolean = false,
            $sortBy: String = "",
            $filterKo: Boolean = false
        ) {
            currentBackground {
                enrich(
                    genes: $genes,
                    filterTerm: $filterTerm,
                    offset: $offset,
                    first: $first,
                    filterFda: $filterFda,
                    sortby: $sortBy,
                    filterKo: $filterKo
                ) {
                    nodes {
                        pvalue
                        adjPvalue
                        oddsRatio
                        nOverlap
                        geneSets {
                            nodes {
                                term
                                id
                                nGeneIds
                            }
                        }
                    }
                }
            }
        }""",
    }

    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json"
    }

    response = requests.post(url, json=query, headers=headers)  # используем json= для автоматической сериализации
    response.raise_for_status()
    res = response.json()

    # Защита от ошибок в структуре ответа (чтоб не выдал None)
    if res.get('data') is None:
        print("Ошибка: сервер вернул data = null. Возможно, ошибка в запросе или нет данных.")
        return pd.DataFrame()

    current_bg = res['data'].get('currentBackground')
    if current_bg is None:
        print("currentBackground is null — возможно, все гены невалидны или список пуст.")
        return pd.DataFrame()

    enrich_result = current_bg.get('enrich')
    if enrich_result is None:
        print("Поле 'enrich' отсутствует в ответе.")
        return pd.DataFrame()

    enrichment_nodes = enrich_result.get('nodes', [])
    if not enrichment_nodes:
        return pd.DataFrame()

    # разворачиваем json в плоский датафрейм
    df = pd.json_normalize(
        enrichment_nodes,
        record_path=['geneSets', 'nodes'],
        meta=['pvalue', 'adjPvalue', 'oddsRatio', 'nOverlap']
    )

    if df.empty:
        return pd.DataFrame()

    # --- Парсим term для извлечения condition annotations ---
    def parse_term(term: str):
        parts = term.split('_')
        direction = "N/A"
        concentration = "N/A"

        # Извлекаем направление (up / down), оно находится после пробела в конце строки
        if ' ' in term:
            main_part, direction = term.rsplit(' ', 1)
        else:
            main_part = term

        # Разбиваем основную часть по '_'
        tokens = main_part.split('_')
        n = len(tokens)

        batch = tokens[0] if n > 0 else "N/A"
        timepoint = tokens[1] if n > 1 else "N/A"
        cellLine = tokens[2] if n > 2 else "N/A"
        batch2 = tokens[3] if n > 3 else "N/A"
        perturbation_raw = tokens[4] if n > 4 else "N/A"

        if n > 5:
            # Концентрация — всё, что после 5-го токена
            concentration = '_'.join(tokens[5:])

        # Очистка perturbation: убираем всё после последнего пробела (если есть направление внутри)
        perturbation_clean = perturbation_raw
        if ' ' in perturbation_raw:
            perturbation_clean = perturbation_raw.split(' ')[0]

        # Обработка CRISPR KO (если строка заканчивается на " KO")
        if perturbation_raw.endswith(" KO"):
            perturbation = perturbation_raw
        else:
            perturbation = perturbation_clean

        return pd.Series({
            'batch': batch,
            'timepoint': timepoint,
            'cellLine': cellLine,
            'batch2': batch2,
            'perturbation': perturbation,
            'concentration': concentration,
            'direction': direction
        })

    # Применяем парсинг
    annotation_cols = df['term'].apply(parse_term)
    df = pd.concat([df, annotation_cols], axis=1)

    # Возвращаем все аннотации + статистику
   # result_cols = [
   #     'term', 'perturbation', 'cellLine', 'timepoint',
   #     'concentration', 'direction',
   #     'pvalue', 'adjPvalue', 'oddsRatio', 'nOverlap'
   # ]
    return df.reset_index(drop=True)


# возможно пригодиться:
# приходится смотреть перекрывание генов из запроса и генов из бекграунда L1000
# (иначе несуществующие в LINCS гены проигнорируются)
'''
def get_l2s2_valid_genes(genes: list):
    query = {
        "operationName": "GenesQuery",
        "variables": {"genes": genes},
        "query": """query GenesQuery($genes: [String]!) {
            geneMap2(genes: $genes) {
                nodes {
                    gene
                    geneInfo {
                        symbol
                    }
                }
            }
        }"""
    }
    response = requests.post(url, json=query)
    response.raise_for_status()
    res = response.json()
    return [g['geneInfo']['symbol'] for g in res['data']['geneMap2']['nodes'] if g.get('geneInfo')]
'''

def get_l2s2_valid_genes(genes: list[str]):
    query = {
    "query": """query GenesQuery($genes: [String]!) {
        geneMap2(genes: $genes) {
            nodes {
                gene
                geneInfo {
                    symbol
                    }
                }
            }
        }""",
    "variables": {"genes": genes},
    "operationName": "GenesQuery"
    }
    
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json"
    }

    response = requests.post(url, data=json.dumps(query), headers=headers)

    response.raise_for_status()
    res = response.json()
    return [g['geneInfo']['symbol'] for g in res['data']['geneMap2']['nodes'] if g['geneInfo'] != None]






#####################
# PARSING ARGUMENTS #
#####################

# Создаем парсер
parser = argparse.ArgumentParser(description="Программа для приветствия пользователя.")
# обязательный аргумент
parser.add_argument("perturbagen_group", type=str, help="The name of the perturbagen group, to which belong observed signatures")
# необязательный флаг/опция
parser.add_argument("-project_path", type=str,
                    default="Users/svetlana/Desktop/CMap_2025-26/",
                    help="Path to the folder with the whole project")
parser.add_argument("-out_folder", type=str,
                    default="L2S2_results_20000",
                    help="The NAME of out-folder, where will be located obtained results of this script. the location of this folder has already been determined by the script.")

# Считываем аргументы
args = parser.parse_args()

pert_grp =  args.perturbagen_group
if args.project_path != "":
    project_folder_path = args.project_folder_path
if args.out_folder != "":
    algorithm_results_folder_name = args.out_folder


directory = f"/{project_folder_path}/{pert_grp}_inhibitors/signatures"
l2s2_res = f"/{project_folder_path}/{pert_grp}_inhibitors/L2S2/{algorithm_results_folder_name}"

with os.scandir(directory) as entries:
    for entry in entries:
        if entry.is_file():
            print(f"Обрабатываю: {entry.name}")
            
            with open(entry.path, 'r', encoding='utf-8') as sgn:
                signature = sgn.read().split("\n")
                signature = [g.strip() for g in signature if g.strip()]  # чистим пробелы и пустые строки
                signature = [g.split("\t")[0] for g in signature]
                
            #print(signature)
            #break
            try:
                print(f"signature length = {len(signature)}")
            
                if signature != []:
                    annotations = get_l2s2_condition_annotations(signature, first=20000)

                    os.chdir(l2s2_res)
                    annotations.to_csv(f"{entry.name[:-4]}.L2S2.tsv", sep='\t', index=False, encoding='utf-8')
                else:
                    print("zero genes input")
                
            except HTTPError as http:
                print(f"error: {http}")
                if http.response.status_code == 413:
                    print(f"signature {entry.name} is too big: {len(signature)} genes")
                continue
            
                
            

            


            
            
            
            
