# %%
# 计算AUC
import numpy as np
from sklearn import metrics
from sklearn.metrics import precision_recall_curve,auc,roc_auc_score
from sklearn.metrics import average_precision_score
import matplotlib.pyplot as plt
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='参数')
parser.add_argument('--input_file', '-i',help='输入原始数据')
parser.add_argument('--input_result', '-r',help='输入预测结果')
parser.add_argument('--output', '-o',help='输出统计结果')
parser.add_argument('--output_result', '-or',help='输出预测结果')

args=parser.parse_args()

input=args.input_file
result=args.input_result
output=args.output
outw=open(args.output_result,'w')


# %%
dict_label={}
with open(input) as f:
    for eachline in f:
        line=eachline.strip().split("\t")
        pos=line[0]+"\t"+line[-1].replace("*","")
        dict_label[pos]=line[1]

# %%
list_true,list_pre=[],[]
#ncorrect, nsamples = 0, 0
with open(result) as f:
    for eachline in f:
        line=eachline.strip().split("\t")
        if "peptide" not in line[0] : 
            pos=line[0]+"\t"+line[2]
            true_label=dict_label[pos]
            pre_label=line[3]
            list_true.append(int(true_label))

            pre=float(line[4].split(',')[-1].replace("]",'').strip())
            list_pre.append(pre)

            w_infor=line[0]+"\t"+line[3]+"\t"+line[2]+"\t"+line[4]+"\t"+true_label
            outw.write("%s\n"%(w_infor))
        elif "peptide" in line[0]:
            w_infor="peptide\tMHC\tHLAtype\tprobability\tReal_Bind"
            outw.write("%s\n"%(w_infor))
outw.close()
# %%
# ====AUC ====
roc_auc = roc_auc_score(list_true,list_pre)
#print("AUC:", roc_auc)

# ====AUC 0.1====
partial_auc = roc_auc_score(list_true,list_pre, max_fpr=0.1)
#print("AUC0.1:", partial_auc)

# ==== Precision-Recall 曲线和 AUC ====
precision, recall, pr_thresholds = precision_recall_curve(list_true,list_pre)
auc_pr = auc(recall, precision)
#print(precision)

# ====输出主要指标====
print("=== 评估结果 ===")
print(f"ROC-AUC (全范围): {roc_auc:.4f}")
print(f"ROC-AUC (FPR ≤ 0.1): {partial_auc:.4f}")
print(f"PR-AUC: {auc_pr:.4f}")
print(f"平均Precision: {np.mean(precision):.4f}")
print(f"平均Recall: {np.mean(recall):.4f}")
w=str(roc_auc)+"\t"+str(partial_auc)+"\t"+str(auc_pr)+"\t"+str(np.mean(recall))+"\t"+str(np.mean(precision))


table="AUC\tAUC0.1\tPR-AUC\tmean-Recall\tmean-precision\n"
out=open(output,'w')
out.write("%s%s\n"%(table,w))
out.close()
# %%



