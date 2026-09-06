import os
import pandas as pd
import csv

### mTOR ###

deg_file_name = "condition_control_case.deseq2.results_filtered.tsv"
de_table = "condition_control_case.deseq2.results_names.tsv"
dir_path = "/data/popov/TDEs"
mTOR_dir = "/data/popov/mTOR_inhibitors/signatures_initial"
mtor_tdes = "/data/popov/mTOR_inhibitors/results_TDEs"

mTOR_dirs = [
    "results_42_59",
    "results_60_79",
    "results_80_99",
    "results_100_119",
    "results_120_135"
    ]

for dr in mTOR_dirs:
    os.chdir(f"{dir_path}/{dr}")
    sgns = os.listdir()
    for sgn in sgns:
        deg_file = os.path.join(dir_path, dr, sgn, deg_file_name)
        degs = pd.read_csv(deg_file, sep="\t")
        tde_file = os.path.join(dir_path, dr, sgn, de_table)
        tde = pd.read_csv(tde_file, sep = "\t")
        
        sgn_out_file = os.path.join(mTOR_dir, f"{sgn}.sgn.txt")
        degs.to_csv(sgn_out_file, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar="\\")
        tde_out_file = os.path.join(mtor_tdes, f"{sgn}.deseq2_table.tsv")
        tde.to_csv(tde_out_file, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar="\\")
    


### PI3K ###

deg_file_name = "condition_control_case.deseq2.results_filtered.tsv"
de_table = "condition_control_case.deseq2.results_names.tsv"
dir_path = "/data/popov/TDEs"
PIEK_dir = "/data/popov/PI3K_inhibitors/signatures_initial"
piek_tdes = "/data/popov/PI3K_inhibitors/results_TDEs"

PIEK_dirs = [
    "results_140_148",
    "results_150_154"
    ]

for dr in PIEK_dirs:
    os.chdir(f"{dir_path}/{dr}")
    sgns = os.listdir()
    for sgn in sgns:
        deg_file = os.path.join(dir_path, dr, sgn, deg_file_name)
        degs = pd.read_csv(deg_file, sep="\t")
        tde_file = os.path.join(dir_path, dr, sgn, de_table)
        tde = pd.read_csv(tde_file, sep = "\t")
        
        out_file = os.path.join(PIEK_dir, f"{sgn}.sgn.txt")
        degs.to_csv(out_file, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar="\\")
        tde_out_file = os.path.join(piek_tdes, f"{sgn}.deseq2_table.tsv")
        tde.to_csv(tde_out_file, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar="\\")


### CDK ###

deg_file_name = "condition_control_case.deseq2.results_filtered.tsv"
de_table = "condition_control_case.deseq2.results_names.tsv"
dir_path = "/data/popov/TDEs"
CDK_dir = "/data/popov/CDK_inhibitors/signatures_initial"
cdk_tdes = "/data/popov/CDK_inhibitors/results_TDEs"

CDK_dirs = [
    "results_180_199",
    "results_200_219",
    "results_220_239",
    "results_240_259"
    ]

for dr in CDK_dirs:
    os.chdir(f"{dir_path}/{dr}")
    sgns = os.listdir()
    for sgn in sgns:
        deg_file = os.path.join(dir_path, dr, sgn, deg_file_name)
        degs = pd.read_csv(deg_file)
        tde_file = os.path.join(dir_path, dr, sgn, de_table)
        tde = pd.read_csv(tde_file, sep = "\t")
        
        out_file = os.path.join(CDK_dir, f"{sgn}.sgn.txt")
        degs.to_csv(out_file, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar="\\")
        tde_out_file = os.path.join(cdk_tdes, f"{sgn}.deseq2_table.tsv")
        tde.to_csv(tde_out_file, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar="\\")

        
