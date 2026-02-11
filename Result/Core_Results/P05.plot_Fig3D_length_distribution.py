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
# -----------------------------
# 颜色映射（按照你指定的颜色）
# -----------------------------
color_map = {
    "NetMHCpan4.1-BA": "#B4D5AF",
    "MixMHCpred2.2": "#91C0C7",
    "mhcflurry2.1.0": "#6E96A7",
    "NetMHCpan4.1-EL": "#6FA961",
    "LKMHCpan": "#D28A39"    
}

# %%
# -----------------------------
# 读取你的 TSV 文件（替换文件名）
# -----------------------------
#df = pd.read_csv("Evaluation_data_set_from_netMHCpan_trim_AUC_peptideLength.xls", sep="\t")
df = pd.read_csv(args.input_file, sep="\t")
# 保证列名一致
# df.columns 应包含："Length", "PR-AUC", "Software"

# 获取所有长度（X 轴顺序）
lengths = sorted(df["Length"].unique())

# 获取软件列表（按你的颜色顺序）
softwares = list(color_map.keys())

# %%
# -----------------------------
# 开始画分组柱状图
# -----------------------------
bar_width = 0.13
cluster_spacing = 0.8
x = np.arange(len(lengths)) * cluster_spacing

plt.figure(figsize=(10, 5))

for i, software in enumerate(softwares):
    sub_df = df[df["Software"] == software] \
                .set_index("Length") \
                .reindex(lengths)          # 按长度排序

    y = sub_df["PR-AUC"].values

    plt.bar(
        x + i * bar_width,
        y,
        width=bar_width,
        label=software,
        color=color_map[software]
    )

# -----------------------------
# 图形美化
# -----------------------------
plt.rcParams["font.family"] = "Arial"
plt.rcParams['font.size'] = 12
plt.xticks(x + bar_width * (len(softwares) / 2), lengths)
plt.xlabel("Peptide Length")
plt.ylabel("PR-AUC")
plt.ylim(0, 1.0)
# plt.legend()
plt.tight_layout()
#plt.savefig("FigC.svg", format="svg")
plt.savefig(args.output_file, format="svg")
plt.show()



# %%



