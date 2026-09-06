#!/usr/bin/env python3

################ L1000_ALGORITHM_RUNNER.py #####################
#
# perturbagen_group  ---  name of signatures perturbagen's group (compulsory option), e.g. MEK, mTOR, PI3K...
# --project_path  ---  path to a directory with the whole project
# --out_folder  ---  NAME of out-folder, where will be located obtained results of this script. the location of this folder has already been determined by the script.
#
################################################################




import json, requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from pprint import pprint
import ast
import sys
import pandas as pd
import os
import time
import argparse


def GSSl1000(file: str, log_file, sgns_dir: str, out_dir: str, headers: json, session: requests.Session):
    print(f"processing file {file}", file=log_file)
    print(f"processing file {file}")

    try:
        sgn = pd.read_csv(sgns_dir + "/" + file, sep="\t")
    except:
        print(f"ERROR: file {file} contains no genes", file=log_file)
        print(f"ERROR: file {file} contains no genes")
        return "empty"
    up_genes = sgn[sgn.iloc[:, 2] > 0].iloc[:, 0].tolist()
    down_genes = sgn[sgn.iloc[:, 2] < 0].iloc[:, 0].tolist()

    print(f"*** N of up-genes: {len(up_genes)}", file=log_file)
    print(f"*** N of down-genes: {len(down_genes)}", file=log_file)
    print(f"*** N of up-genes: {len(up_genes)}")
    print(f"*** N of down-genes: {len(down_genes)}")

    payload = {
        'up_genes': up_genes,
        'down_genes': down_genes
    }

    #  try:
    # response = requests.post(L1000FWD_URL + 'sig_search', json=payload)
    with session.post(L1000FWD_URL + 'sig_search', headers=headers, json=payload, stream=True) as response:
        if response.status_code != 200:
            results = pd.DataFrame()
            results.to_csv(out_dir + "/" + file[:-4] + "l1000.tsv", index=False, sep="\t")
            print(f"POST: status code IS NOT 200: {response.status_code}", file=log_file)
            print(f"POST: status code IS NOT 200: {response.status_code}")
            return "not 200"

        time.sleep(2)

        result_ids = list(response.json().values())
        for result_id in result_ids:
            #print("\n")
            #print(f"RESULT_ID: {result_id}", file=log_file)
            #print(f"RESULT_ID: {result_id}")
            #print("\n")
            try:
                # response = requests.get(f"https://maayanlab.cloud/l1000fwd/result/graph/{result_id}?n=100") # !!!!!!!!!!!!!!!!
                res_url = f"https://maayanlab.cloud/l1000fwd/result/graph/{result_id}?n=100"
                with session.get(res_url, headers=headers, timeout=30, stream=True) as response:
                    if response.status_code != 200:
                        results = pd.DataFrame()
                        results.to_csv(out_dir + "/" + file[:-4] + "l1000.tsv", index=False, sep="\t")
                        print(f"GET: status code IS NOT 200: {response.status_code}", file=log_file)
                        print(f"GET: status code IS NOT 200: {response.status_code}")
                        continue

                    print(f"Size of l1000 results: {sys.getsizeof(response.content)}", file=log_file)
                    print(f"Size of l1000 results: {sys.getsizeof(response.content)}")

                    raw_data = response.content
                    results = pd.DataFrame(json.loads(raw_data))

                    print(f"Number of findings: {len(results)}", file=log_file)
                    print(f"Number of findings: {len(results)}")

                    results.to_csv(out_dir + "/" + file[:-4] + "l1000.tsv", index=False, sep="\t")

            except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
                print(f"ERROR downloading result {result_id}: {e}", file=log_file)
                print(f"ERROR downloading result {result_id}: {e}")
                continue

    return "ok"


#####################
# PARSING ARGUMENTS #
#####################

# Создаем парсер
parser = argparse.ArgumentParser(description="Программа для приветствия пользователя.")
# обязательный аргумент
parser.add_argument("perturbagen_group", type=str,
                    help="The name of the perturbagen group, to which belong observed signatures")
# необязательный флаг/опция
parser.add_argument("-project_path", type=str,
                    default="Users/svetlana/Desktop/CMap_2025-26/",
                    help="Path to the folder with the whole project")
parser.add_argument("-out_folder", type=str,
                    default="L1000_results_with_empty_files",
                    help="The NAME of out-folder, where will be located obtained results of this script. the location of this folder has already been determined by the script.")

# Считываем аргументы
args = parser.parse_args()

pert_grp = args.perturbagen_group
if args.project_path != "":
    project_folder_path = args.project_folder_path
if args.out_folder != "":
    algorithm_results_folder_name = args.out_folder



L1000FWD_URL = 'https://maayanlab.cloud/l1000fwd/'


sgns_dir = f"/{project_folder_path}/{pert_grp}_inhibitors/signatures"
out_dir = f"/{project_folder_path}/{pert_grp}_inhibitors/L1000/{algorithm_results_folder_name}"
logs = f"/{project_folder_path}/{pert_grp}_inhibitors/L1000/l1000.logs.txt"
txt_files = [f for f in os.listdir(sgns_dir) if f.endswith('.txt')]
sgns = sorted(txt_files, key=lambda x: int(x.split(".")[0]))

max_retries = 5

with open(logs, "a") as log_file:

   # sgns = sorted(os.listdir(sgns_dir))
   # sgns = sgns[34:35]

    # 1. Используем сессию с автоматическими повторами
    print("\n\n", file=log_file)
    print(f"Creating new requests Session (max_retries={max_retries})", file=log_file)
    print(f"Creating new requests Session (max_retries={max_retries})")
    session = requests.Session()
    retry_strategy = Retry(
        total=max_retries,
        backoff_factor=2,  # Время ожидания между повторами (2, 4, 8... секунд)
        raise_on_status=False
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount('https://', adapter)
    headers = {'User-Agent': 'Mozilla/5.0'}

    i = 0
    while i < len(sgns):
        file = sgns[i]
        try:
            status = GSSl1000(file, log_file, sgns_dir, out_dir, headers, session)
            if status == "empty" or status == "not 200":
                i+=1
                continue
            elif status == "ok":
                print(f"SUCCESS: {sgns[i]}", file=log_file)
                print(f"SUCCESS: {sgns[i]}")
                i += 1
                continue

        except requests.exceptions.ConnectionError as ce:
            print(f"CONNECTION ERROR ({ce}): TRYING AGAIN: {file}", file=log_file)
            print(f"CONNECTION ERROR ({ce}): TRYING AGAIN: {file}")
            continue # trying sgns[i] again




