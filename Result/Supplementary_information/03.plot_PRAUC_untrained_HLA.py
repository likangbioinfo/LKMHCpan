# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='参数')
parser.add_argument('--input_file', '-i', help='输入文件')
parser.add_argument('--output_file', '-o', help='输出画图文件')
args=parser.parse_args()

# %%
plt.figure(figsize=(6, 4))
plt.rcParams["font.family"] = "Arial"
plt.rcParams['font.size'] = 12

# %%
# -----------------------------
# 颜色映射（按照你指定的颜色）
# -----------------------------
color_map = {
    "LKMHCpan": "#D28A39"    
}

# %%
# ---------------------------------------------------------
# 读取输入文件（TSV）
# ---------------------------------------------------------
#df = pd.read_csv("Evaluation_data_set_from_netMHCpan_trim_AUC_HLA.xls", sep="\t")
df = pd.read_csv(args.input_file, sep="\t")
df_sort=df.sort_values(by='PR-AUC', ascending=False)
plt.bar(df_sort['Length'], df_sort['PR-AUC'],color='#D28A39',width=0.6)
plt.xticks(rotation=45, ha='right')   # x 轴

plt.ylabel("PR-AUC")
plt.ylim(0, 1.0)

#plt.savefig("FigD.svg", format="svg")
plt.savefig(args.output_file, format="svg")
plt.show()


