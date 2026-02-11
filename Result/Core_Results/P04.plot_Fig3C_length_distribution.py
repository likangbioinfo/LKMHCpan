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
plt.rcParams["font.family"] = "Arial"
plt.rcParams['font.size'] = 12

# %%
# -----------------------------
# 颜色映射（按照你指定的颜色）
# -----------------------------
color_map = {
    "NetMHCpan4.1-BA": "#B4D5AF",
    "MixMHCpred2.1": "#91C0C7",
    "mhcflurry2.0.2": "#6E96A7",
    "NetMHCpan4.1-EL": "#6FA961",
    "LKMHCpan": "#D28A39"    
}

# %%
# ---------------------------------------------------------
# 读取输入文件（TSV）
# ---------------------------------------------------------
#df = pd.read_csv("Evaluation_data_set_from_netMHCpan_trim_AUC_HLA.xls", sep="\t")
df = pd.read_csv(args.input_file, sep="\t")
# 唯一的 HLA 类型（x 轴）
hlas = df["HLA"].unique().tolist()

# 软件列表按你定义顺序
softwares = list(color_map.keys())
print(softwares)


# %%
# ---------------------------------------------------------
# 开始绘图
# ---------------------------------------------------------
num_hla = len(hlas)
num_sw = len(softwares)
bar_width = 0.8 / num_sw
x = np.arange(num_hla)

plt.figure(figsize=(26, 4))

# 为每个软件绘制柱子
for i, sw in enumerate(softwares):
    sub_df = df[df["Software"] == sw].set_index("HLA")
    y = sub_df.reindex(hlas)["PR-AUC"].values  # 按 HLA 顺序提取
    
    plt.bar(
        x + i * bar_width,
        y,
        width=bar_width,
        label=sw,
        color=color_map[sw]
    )
# ---------------------------------------------------------
# 图形格式美化
# ---------------------------------------------------------
plt.xticks(x + bar_width * num_sw / 2, hlas, rotation=45)
plt.ylabel("PR-AUC")
plt.ylim(0, 1.0)

plt.xlim(-bar_width-0.2, num_hla - 1 + num_sw * bar_width+0.2)
# plt.legend(ncol=5)
plt.tight_layout()

#plt.savefig("FigD.svg", format="svg")
plt.savefig(args.output_file, format="svg")
plt.show()


