import os
import glob
from platform import python_version
import sys

import pandas as pd
import numpy as np
from pandas.core.interchange.dataframe_protocol import DataFrame
import matplotlib.pyplot as plt


inh_grp = "mTOR"#sys.argv[1]

lincs = pd.read_csv("/Users/svetlana/Desktop/CMap_2025-26/GSE92742_Broad_LINCS_gene_info.txt", sep="\t")

try:
    os.chdir(f"/Users/svetlana/Desktop/CMap_2025-26/{inh_grp}_inhibitors/signatures")
except:
    print(f"No signatures` files for {inh_grp} inhibitors")

sgns = glob.glob("*.txt")
sgns = sorted(sgns, key=lambda x: int(x.split(".")[0]))

lm = set(lincs[lincs["pr_is_lm"] == 1]["pr_gene_symbol"])
print(f"N of LM: {len(lm)}")
all_lincs = set(lincs["pr_gene_symbol"])
print(f"N of all lincs genes: {len(all_lincs)}")


all_genes = set()
portions_per_sgn = {
    "sgn_lm": [],
    "sgn_bing": [],
    "sgn_notlincs": []
}
for file in sgns:
    print(f"PROCESSING {file}")
    try:
        sgn = pd.read_csv(file, sep="\t")
    except:
        print(f"file {file} is empty")
        continue
    sgn = set(sgn.iloc[:,0])
    for j in sgn:
        all_genes.add(j)
    sgn_lm = len(sgn & lm)
    sgn_bing = len(sgn & all_lincs) - sgn_lm
    sgn_notlincs = len(sgn) - sgn_bing - sgn_lm
    portions_per_sgn["sgn_lm"].append(sgn_lm)
    portions_per_sgn["sgn_bing"].append(sgn_bing)
    portions_per_sgn["sgn_notlincs"].append(sgn_notlincs)


        
print(f"N of all genes: {len(all_genes)}")


print([len(i) for i in portions_per_sgn.values()])


# is LM, is in set
TT = lm & all_genes
print(f"N of TT (is LM, is in set): {len(TT)}")

# not LM, is in set
FT = all_genes - TT
print(f"N of FT (not LM, is in set): {len(FT)}")

# is LM, not in set
TF = lm - TT
print(f"N of TF (is LM, not in set): {len(TF)}")


## what's in lincs
# in lincs & in set
TT_lincs = all_lincs & all_genes
print(f"N of genes, which are in LINCS: {len(TT_lincs)}")

# in lincs & not in set
TF_lincs = all_lincs - TT_lincs
print(f"N of set genes, which are not in LINCS: {len(TF_lincs)}")

# not in lincs & in set
FT_lincs = all_genes - TT_lincs
print(f"N of LINCS genes, which are not in set: {len(FT_lincs)}")
    
    
    
#print(portions_per_sgn)




portions_per_sgn = pd.DataFrame(portions_per_sgn)
normalized_portions = portions_per_sgn.div(portions_per_sgn.sum(axis=1), axis=0)
sgns.sort(key=lambda x: int(x.split(".")[0]))
normalized_portions.index = sgns

# check
# print(normalized_portions.sum(axis=1))

df = normalized_portions

#ax = df.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='viridis')

# 3. Настройка оси Y в проценты (от 0 до 100)
# Мы меняем метки оси, чтобы они отображались как "%"
#yticks = ax.get_yticks()

#plt.close() # Закроем предыдущий график, чтобы сделать это правильно через нормализацию данных

# --- ПРАВИЛЬНЫЙ ПОДХОД ЧЕРЕЗ НОРМАЛИЗАЦИЮ ДАННЫХ ---

df_norm = df
df_norm.index = [f.removesuffix(".sgn.txt") for f in df_norm.index]


fig, ax = plt.subplots(figsize=(10, 6))
df_norm.plot(kind='bar', stacked=True, ax=ax, colormap='plasma', width=0.8)
ax.tick_params(axis='x', labelsize=7)

# Форматируем ось Y в проценты
from matplotlib.ticker import PercentFormatter
ax.yaxis.set_major_formatter(PercentFormatter(1.0)) # 1.0 означает 100%

ax.set_ylabel('Доля (%)', fontsize=12)
ax.set_title('Состав сигнатур: какие гены есть в LINCS, каких нет', fontsize=14)
ax.legend(title='подмножества', loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylim(0, 1.05) # Чуть выше 100% для красоты
ax.grid(axis='y', linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()
try:
    os.chdir(f"/Users/svetlana/Desktop/CMap_2025-26/{inh_grp}_inhibitors/visualization")
except:
    print("You forgot to create a directory for plots")
plt.savefig('sgns&LINCS.png')


