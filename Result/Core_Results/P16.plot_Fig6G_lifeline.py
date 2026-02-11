# %%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

###根据具有免疫原性肽的阈值，计算neoantigen load
###使用neoantigen load对PD vs SD/CR/PR 绘制箱线图，找到阈值,进行生存分析
parser = argparse.ArgumentParser(description='参数')
parser.add_argument('--input_neo', '-in', help='输入样本的neoantigen load')
parser.add_argument('--input_os', '-io', help='输入样本的生存期')
parser.add_argument('--output', '-o', help='输出分组结果文件')
parser.add_argument('--outputplot', '-op', help='输出绘图结果文件')
args=parser.parse_args()

# %%
raw_file=args.input_neo
clinic_file=args.input_os
prefix=(args.outputplot).replace(".svg","")
# %%
raw_df=pd.read_csv(raw_file,sep="\t")
dict_tyep=raw_df[raw_df['Type']=='Lung'].set_index('Sample')['Neoantigen'].to_dict()

# %%
clinic_df=pd.read_csv(clinic_file,sep="\t")

clinic_df_fit=clinic_df[clinic_df['sample'].isin(dict_tyep.keys())]

# %%
clinic_df_fit['neoantigen_load'] = clinic_df_fit['sample'].map(dict_tyep)
list_cd8_pos=clinic_df_fit[clinic_df_fit['Harmonized_Confirmed_BOR']!='PD']['neoantigen_load'].tolist()
list_cd8_neg=clinic_df_fit[clinic_df_fit['Harmonized_Confirmed_BOR']=='PD']['neoantigen_load'].tolist()
print('R:',len(list_cd8_pos))
print('NR:',len(list_cd8_neg))

# %%
#绘制新抗原负荷分布图（条形图）：整体新抗原负荷
prefix_load=prefix+"_neoload.svg"
plt.figure(figsize=(3, 4))
plt.hist(clinic_df_fit['neoantigen_load'], bins=10, edgecolor='black')  
plt.xlabel('neoantigen load')
plt.ylabel('Number')
plt.title('The Distribution of neoantigen load (Overall)')
plt.savefig(prefix_load)
plt.show()

#绘制新抗原负荷分布图（条形图）：Response样本新抗原负荷
clinic_df_fit_response=clinic_df_fit[clinic_df_fit['Harmonized_Confirmed_BOR']!='PD']   
prefix_load_response=prefix+"_neoload_response.svg"
plt.figure(figsize=(3, 4))
plt.hist(clinic_df_fit_response['neoantigen_load'], bins=10, edgecolor  ='black')  
plt.xlabel('neoantigen load')
plt.ylabel('Number')
plt.title('The Distribution of neoantigen load (Response)')
plt.savefig(prefix_load_response)
plt.show()

# %%
#绘制绘制新抗原负荷箱线图
prefix_Response=prefix+"_Response.svg"
data = [list_cd8_pos, list_cd8_neg]
plt.figure(figsize=(2, 4))
bp = plt.boxplot(data, patch_artist=True,widths=0.4)  
# 计算均值和中位数
mean_list1 = np.mean(list_cd8_pos)
median_list1 = np.median(list_cd8_pos)
mean_list2 = np.mean(list_cd8_neg)
median_list2 = np.median(list_cd8_neg)

# 为箱线图添加颜色
colors = ['lightblue', 'lightgreen']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# 获取 x 轴刻度位置
positions = np.arange(1, len(data) + 1)

# 标注均值和中位数
for i, (mean, median) in enumerate([(mean_list1, median_list1), (mean_list2, median_list2)]):
    # 标注均值
    plt.text(positions[i]+0.4, mean, f'Mean: {mean:.2f}', ha='center', va='bottom')
    # 标注中位数
    plt.text(positions[i]+0.4, median, f'Median: {median:.2f}', ha='center', va='top')


# 添加标题和标签
#plt.title('Grouped Boxplot with Mean and Median')
plt.xlabel('Response Type')
plt.ylabel('Neoantigen load')
plt.xticks(positions, ['Response', 'Non-Response'])  
plt.savefig(prefix_Response)

# 显示网格
#plt.grid(True)
plt.show()

# %%
# 定义条件列表
""" conditions=[(clinic_df_fit['neoantigen_load']>0 ) & (clinic_df_fit['neoantigen_load']<=200),
           (clinic_df_fit['neoantigen_load']>200) & (clinic_df_fit['neoantigen_load']<=260),
           (clinic_df_fit['neoantigen_load']>260)
            ] """
medain_cut=np.median(clinic_df_fit['neoantigen_load'].tolist()) ##使用中位数
#medain_cut=200 ##固定一个常数用于新抗原负荷分组High/Low
#medain_cut=np.percentile(np.array(clinic_df_fit['neoantigen_load'].tolist()),75) ##使用上四分位数进行分组
print('cut off:',medain_cut)
conditions=[(clinic_df_fit['neoantigen_load']>0 ) & (clinic_df_fit['neoantigen_load']<=medain_cut),
           (clinic_df_fit['neoantigen_load']>medain_cut)
            ]
# 定义每个条件对应的值
choices = ['Group.L', 'Group.H']

# 使用np.select()创建新列
clinic_df_fit['group'] = np.select(conditions, choices, default='other')  # 'other' 用于不满足任何条件的情况
clinic_df_fit_clean = clinic_df_fit.dropna(subset=['Harmonized_OS_Event', 'Harmonized_OS_Days'])
#clinic_df_fit_clean.to_csv('Table_mutant_merges_peptide_all_predict_8_11_aff_0.5_group.tsv',sep="\t",index=False)
clinic_df_fit_clean.to_csv(args.output,sep="\t",index=False)

# %%
cmd_exe="python3.11 P17.lifeline_neoantige.py -i "+args.output+" -o" +args.outputplot
#print(cmd_exe)
os.system(cmd_exe)
print('#########Finish#################')


