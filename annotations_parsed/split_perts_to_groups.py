from json.decoder import NaN

import pandas as pd
import os
import json
import numpy as np

experiments = pd.read_csv("/Users/svetlana/Desktop/CMap_2025-26/Методы поиск лекарств - Сигнатуры.tsv", sep="\t")
perts = pd.read_csv("/Users/svetlana/Desktop/CMap_2025-26/Методы поиск лекарств - Лекарства по механизму rna-seq.tsv", sep="\t")

pert_ann = perts.iloc[:,[0,1]][0:41]
pert_ann = pert_ann.rename(columns={pert_ann.columns[0]: 'pert', "Тип": "group"})

os.chdir("/Users/svetlana/Desktop/CMap_2025-26")
pert_ann.to_csv("perturbagens_annotation.csv", index=False)
# renamed "perturbagens_annotation.csv" to "perturbagens_MoA_groups.csv"
# and manually added some perts (synonyms from the brackets)
pert_ann = pd.read_csv("perturbagens_MoA_groups.csv")

group = set(pert_ann["group"])
all_perts = set(pert_ann["pert"])

experiments["Signature"]
experiments.columns
all_drugs = set(experiments["perturbation"])
ctrls = all_drugs - all_perts
ctrls

# fixing mistakes of the google table
ctrls.discard("selumetinib")
ctrls.discard("Refamatinib")
ctrls

file_sgns = "/Users/svetlana/Desktop/CMap_2025-26/Методы поиск лекарств - Сигнатуры.tsv"
sgn_pert_match = {}
with open(file_sgns, "r") as f:
    strings = f.readlines()
    strings = strings[1:]
    n_sgn = 0
    pert = ""
    for line in strings:
        line = line.rstrip("\n").split("\t")
        if line[1] != '':
            n_sgn = line[1]
            if line[5] not in ctrls:
                pert = line[5]
                sgn_pert_match[n_sgn] = pert
            else:
                continue
        else:
            if line[5] not in ctrls:
                pert = line[5]
                sgn_pert_match[n_sgn] = pert
            else:
                continue

sgn_pert_match
len(sgn_pert_match)

# fixing google table mistakes
sgn_pert_match["22"] = "Selumetinib"
sgn_pert_match["104"] # = Rapamycin
sgn_pert_match["9"] = "Refametinib" # only this sgn is with Refametinib, so i changed Refamatinib -> Refametinib

# checking for mismatches in the order of sgns in initial google table and generated with the cycle list
len(sgn_pert_match) # 291
sgns_nums = pd.read_csv(file_sgns, sep="\t").iloc[:,1].dropna()
sgns_nums = np.array([int(i) for i in sgns_nums])
len(sgns_nums) # 291
(np.array([int(i) for i in sgn_pert_match.keys()]) == sgns_nums).all() # True

# checking, which sgn is absent
np.argmax(((np.array([int(i) for i in sgn_pert_match.keys()]) == np.array(list(range(1,292)))) - 1).astype(bool))
# = 148, thus there is no sgn number 149

df = pd.DataFrame({"sgn_number": list(sgn_pert_match.keys()), "perturbation": list(sgn_pert_match.values())})
os.chdir("/Users/svetlana/Desktop/CMap_2025-26")
df.to_csv("sgn-perturbagen_annotation.tsv", sep="\t", index=False)

### annotation is ready!

# now lets construct a json to easily obtain all sgns numbers for each pert group
pert_group_match = pert_ann.set_index("pert")["group"].to_dict()
pert_group_match["Refametinib"] = "MEK inhibitors"
group_sgn_match = {}
with open("sgn-perturbagen_annotation.tsv", "r") as f:
    f.readline()
    strings = f.readlines()
    for line in strings:
        line = line.rstrip("\n").split("\t")
        if pert_group_match[line[1]] not in group_sgn_match.keys():
            group_sgn_match[pert_group_match[line[1]]] = [line[0]]
        else:
            group_sgn_match[pert_group_match[line[1]]].append(line[0])

group_sgn_match.keys()

with open('group-sgn_annotation.json', 'w', encoding='utf-8') as f:
    json.dump(group_sgn_match, f, ensure_ascii=False, indent=4)

with open('group-sgn_annotation.json', 'r', encoding='utf-8') as f:
    group_sgn_match = json.load(f)

group_sgn_match["MEK inhibitors"] # 1 ---- 41
group_sgn_match["mTOR inhibitors"] # 42 ---- 135
group_sgn_match["PI3K inhibitors"] # 136 --(-149)-- 170
group_sgn_match["CDK inhibitors"] # 171 ---- 292

