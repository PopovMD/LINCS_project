

############################# SYNONYMS_FINDER.py ###################################
#                                                                                  #
# --project  ---  option to determine the project's folder name                    #
# --inhibitor  ---  option to determine the class of inhibitors (its targets)      #
#                                                                                  #
####################################################################################


import requests
from urllib.parse import quote

def get_pubchem_synonyms(compound_name):
    """
    Получает список синонимов для вещества по его названию из PubChem.
    """
    # Кодируем название для безопасной передачи в URL (на случай пробелов или спецсимволов)
    encoded_name = quote(compound_name)
    
    # Формируем URL запроса
    # Используем endpoint: /compound/name/{name}/synonyms/JSON
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/synonyms/JSON"
    
    try:
        response = requests.get(url, timeout=10)
        
        # Проверка статуса ответа
        if response.status_code == 200:
            data = response.json()
            
            # Структура ответа PubChem:
            # {'InformationList': {'Information': [{'CID': ..., 'Synonym': [...]}]}}
            info_list = data.get('InformationList', {}).get('Information', [])
            
            if not info_list:
                print(f"Информация для '{compound_name}' не найдена.")
                return []
        
            synonyms = info_list[0].get('Synonym', [])
            print("FOUND SOME:")
            return synonyms
            
        elif response.status_code == 404:
            print(f"Вещество '{compound_name}' не найдено в PubChem (Ошибка 404).")
            return []
        else:
            print(f"Ошибка запроса: Статус {response.status_code}")
            return []
            
    except requests.exceptions.RequestException as e:
        print(f"Произошла ошибка при соединении с API: {e}")
        return []



import json
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--project", type=str, default="/Users/svetlana/Desktop/CMap_2025-26/")
parser.add_argument("--inhibitor", type=str, default="MEK")
args = parser.parse_args()


perts_file = str(args.project) + "Signatures_" + str(args.inhibitor) + "-inhibitors/" + str(args.inhibitor) + "_perturbagens_uniq.tsv" #"/Users/svetlana/Desktop/CMap_2025-26/Signatures_MEK-inhibitors/MEK_perturbagens_uniq.tsv"
pert_syn = str(args.project) + "Signatures_" + str(args.inhibitor) + "-inhibitors/" + str(args.inhibitor) + "_perturbagens_synonyms.json"

dic = {}

with open(perts_file, "r") as prt:
    p = prt.readlines()
    for line in p:
        synonyms = get_pubchem_synonyms(line.strip())
        dic[line.strip()] = synonyms

with open(pert_syn, "w", encoding="utf-8") as file:
    json.dump(dic, file, ensure_ascii=False, indent=4)        
        


