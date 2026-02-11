# %%
#######################################
##      File Name: 06.stat_HLA_combination_plot.py
##      Author: liqing
##      Mail: liqing@likanglife.com 
##      Created Time: 2023年04月19日 星期三 17时28分19秒
######################################################
#!/usr/bin/python
# -*- coding:utf-8 -*-

#统计HLA&不同肽段长度组合的对应的数量，并画热图

# %%
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

# %%

parser = argparse.ArgumentParser(description='参数')
parser.add_argument('--hlafile', '-a', help='输入hla伪序列文件')
parser.add_argument('--evalfile', '-e', help='输入Evaluation数据集')
parser.add_argument('--inputfile', '-i', help='输入用于训练的数据集')
args=parser.parse_args()


hlafile=args.hlafile
evalfile=args.evalfile
inputfile=args.inputfile

# %% [markdown]

# %%
dict_hla={}
with open(hlafile) as f:
    for eachline in f:
        line=eachline.split(",")
        if "normalized_allele" not in eachline:
            dict_hla[line[0]]=line[0]

# %%
dict_combine={}
dict_hla_number={}
with open(inputfile) as f:
    for eachline in f:
        if "peptide" not in eachline:
            line=eachline.split(",")
            pep,lab,hla=line[:]
            hla=hla.strip()
            pep_len=len(pep)
            if lab=="1":
                if hla in dict_hla:
                    infor=hla+"-"+str(pep_len)

                    if hla not in dict_hla_number:
                        dict_hla_number[hla]=1
                    else:
                        dict_hla_number[hla]+=1

                    if infor not in dict_combine:
                        dict_combine[infor]=1
                    else:
                        dict_combine[infor]+=1


# %%
x_len,y_len=[],[]
for q in dict_combine:
    if int(dict_combine[q])>=100:
        hlatype=q.replace("HLA-","")
        x_len.append(hlatype+"\n"+str(dict_combine[q]))
        y_len.append(int(dict_combine[q]))
# 将列表转换为数组并更改形状
arr_len_x = np.array(x_len).reshape((21, 20))
arr_len_y = np.array(y_len).reshape((21, 20))       

# %%
len(y_len)

# %% [markdown]
# HLA&肽段长度 类型及其对应的数量

# %%
np.array(y_len).shape

# %%
x,y=[],[]
for m in dict_hla_number:
    #print(m,dict_hla_number[m])
    hlatype=m.replace("HLA-","")
    x.append(hlatype+"\n"+str(dict_hla_number[m]))
    y.append(int(dict_hla_number[m]))
# 将列表转换为数组并更改形状
arr_x = np.array(x).reshape((7, 17))
arr_y = np.array(y).reshape((7, 17))

# %%
len(x)

# %% [markdown]
# 画图：HLA类型及其对应的数量

# %%
##画图: HLA亚型和其对应的数量


# 绘制热图
fig = plt.figure(figsize=(10, 10))

#plt.imshow(arr_y, cmap='hot', interpolation='nearest')
plt.imshow(arr_y,  interpolation='nearest')


# 在热图上标注信息
for i in range(arr_y.shape[0]):
    for j in range(arr_y.shape[1]):
        plt.text(j, i, arr_x[i, j],
                 ha="center", va="center", color="black",fontsize=6)

# 添加色度条，并设置其大小和长度
cbar = plt.colorbar(shrink=0.3)
#cbar.ax.tick_params(labelsize=12, length=5, width=2)

# 设置坐标轴标签和标题
plt.ylabel('')
plt.xlabel('Training Dataset',fontsize=12)
plt.title('The heatmap for positive Training Dataset',fontsize=14)
fig.savefig("Training_dataset_HLAtype_Number.svg")
#outplot=args.output_file
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

# 设置图像的DPI值
fig.set_dpi(600)  # 设置为300 DPI
#plt.savefig("")
# 显示图像
plt.show()

# %% [markdown]
# 画图：HLA&肽段长度 类型及其对应的数量

# %%
##画图: HLA亚型和其对应的数量

# 绘制热图
fig = plt.figure(figsize=(20, 27))

#plt.imshow(arr_y, cmap='hot', interpolation='nearest')
plt.imshow(arr_len_y,  interpolation='nearest')


# 在热图上标注信息
for i in range(arr_len_y.shape[0]):
    for j in range(arr_len_y.shape[1]):
        plt.text(j, i, arr_len_x[i, j],
                 ha="center", va="center", color="black",fontsize=8)

# 添加色度条，并设置其大小和长度
cbar = plt.colorbar(shrink=0.3)
#cbar.ax.tick_params(labelsize=12, length=5, width=2)

# 设置坐标轴标签和标题
plt.ylabel('')
plt.xlabel('Training Dataset',fontsize=14)
plt.title('The heatmap for positive Training Dataset',fontsize=20)
fig.savefig("Training_dataset_HLA_peptidesLength_Number.svg")
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# 设置图像的DPI值
fig.set_dpi(600)  # 设置为300 DPI
# 显示图像
plt.show()

# %% [markdown]
# * 统计 Evaluation Dataset中的HLA、HLA&肽段长度

# %%
dict_combine_eval={}
dict_hla_number_eval={}
with open(evalfile) as f:
    for eachline in f:
        if "peptide" not in eachline:
            line=eachline.split()
            pep,lab,hlaseq,hla=line[:]
            pep_len=len(pep)
            if lab=="1":
                if hla in dict_hla:
                    infor=hla+"-"+str(pep_len)

                    if hla not in dict_hla_number_eval:
                        dict_hla_number_eval[hla]=1
                    else:
                        dict_hla_number_eval[hla]+=1

                    if infor not in dict_combine_eval:
                        dict_combine_eval[infor]=1
                    else:
                        dict_combine_eval[infor]+=1

# %%
print(len(dict_combine_eval))
print(len(dict_hla_number_eval))

# %%
x_len_eval,y_len_eval=[],[]
for q in dict_combine_eval:
    if int(dict_combine_eval[q])>=0:
        hlatype=q.replace("HLA-","")
        x_len_eval.append(hlatype+"\n"+str(dict_combine_eval[q]))
        y_len_eval.append(int(dict_combine_eval[q]))
# 将列表转换为数组并更改形状
arr_len_x_eval = np.array(x_len_eval).reshape((7, 35))
arr_len_y_eval= np.array(y_len_eval).reshape((7, 35))       

# %%
x_eval,y_eval=[],[]
for m in dict_hla_number_eval:
    #print(m,dict_hla_number[m])
    hlatype=m.replace("HLA-","")
    x_eval.append(hlatype+"\n"+str(dict_hla_number_eval[m]))
    y_eval.append(int(dict_hla_number_eval[m]))
# 将列表转换为数组并更改形状
arr_x_eval = np.array(x_eval).reshape((7, 5))
arr_y_eval = np.array(y_eval).reshape((7, 5))

# %%
##画图: HLA亚型和其对应的数量

# 绘制热图
fig = plt.figure(figsize=(6, 3))

#plt.imshow(arr_y, cmap='hot', interpolation='nearest')
plt.imshow(arr_y_eval,  interpolation='nearest')


# 在热图上标注信息
for i in range(arr_y_eval.shape[0]):
    for j in range(arr_y_eval.shape[1]):
        plt.text(j, i, arr_x_eval[i, j],
                 ha="center", va="center", color="black",fontsize=4)

# 添加色度条，并设置其大小和长度
cbar = plt.colorbar(shrink=0.3)
#cbar.ax.tick_params(labelsize=12, length=5, width=2)

# 设置坐标轴标签和标题
plt.ylabel('')
plt.xlabel('Evaluation Dataset',fontsize=8)
plt.title('The heatmap for positive Evaluation Dataset',fontsize=8)
fig.savefig("Evaluation_dataset_HLAtype_Number.svg")
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)

# 设置图像的DPI值
fig.set_dpi(600)  # 设置为300 DPI
# 显示图像
plt.show()

# %%
##画图: HLA亚型和其对应的数量

# 绘制热图
fig = plt.figure(figsize=(26, 12))

#plt.imshow(arr_y, cmap='hot', interpolation='nearest')
plt.imshow(arr_len_y_eval,  interpolation='nearest')


# 在热图上标注信息
for i in range(arr_len_y_eval.shape[0]):
    for j in range(arr_len_y_eval.shape[1]):
        plt.text(j, i, arr_len_x_eval[i, j],
                 ha="center", va="center", color="black",fontsize=7)

# 添加色度条，并设置其大小和长度
cbar = plt.colorbar(shrink=0.15,location="right")
#cbar.ax.tick_params(labelsize=12, length=5, width=2)

# 设置坐标轴标签和标题
plt.ylabel('')
plt.xlabel('Evaluation Dataset',fontsize=14)
plt.title('The heatmap for positive Evaluation Dataset',fontsize=20)
fig.savefig("Evaluation_dataset_HLA_peptidesLength_Number.svg")
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# 设置图像的DPI值
fig.set_dpi(600)  # 设置为300 DPI
# 显示图像
plt.show()


