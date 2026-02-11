# %%
import os
import sys
import pandas as pd
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
import warnings
import argparse
warnings.filterwarnings("ignore")  
#肺癌+黑色素瘤
#统计每个样本的新抗原数量

# %%
parser = argparse.ArgumentParser(description='参数')
parser.add_argument('--input_file', '-i', help='输入原始的结果')
parser.add_argument('--input_cutoff', '-c', help='输入阈值')
parser.add_argument('--output', '-o', help='输出绘图结果文件')
parser.add_argument('--intput_type', '-t', help='输出软件类型')
parser.add_argument('--output_result', '-or', help='输出neoantigen load结果文件')
args=parser.parse_args()
####读原始结果文件
df_raw_file=args.input_file
# %%
####读原始结果文件
cutoff=float(args.input_cutoff)
df_data_raw_all=pd.read_csv(df_raw_file,sep="\t")
intput_type=args.intput_type
#df_data_raw=df_data_raw_all[df_data_raw_all['LKMHCpan_pre']>=cutoff]
df_data_raw=df_data_raw_all[df_data_raw_all[intput_type]>=cutoff]
df_count_dict=df_data_raw['ID'].value_counts().to_dict()
df_response_dict=df_data_raw.set_index('ID')['Immunity'].to_dict()

# %%
print('############Start#################')
print(args.input_file)
print('阈值:',cutoff)
print('检出新抗原的样本总数:',len(df_response_dict))

# %%
list_sample=[]
for sample in df_response_dict:
    response=df_response_dict[sample]
    number=df_count_dict[sample]
    #print(sample,number,response)
    if sample.startswith('SU'):
        #print(sample)
        type="Lung"
    else:
        type="Melanoma"
    infor=type+"\t"+str(sample)+"\t"+str(number)+"\t"+response
    list_sample.append(infor)
#print(list_sample)
# %%
df_result=pd.DataFrame([x.split('\t') for x in list_sample])
df_result.columns=['Type','Sample','Neoantigen','Response']
df_result.to_csv(args.output_result,sep="\t",index=False)
#print(df_result.head())
# %% [markdown]
# 绘制黑色素瘤的新抗原和响应的关系图

# %%
###anti PD-1 experienced
str_mela="3754,3909,3919,3936,4044,4084,4085,4105,4120,4135,4142,4147,4154,4165,4169,4178,4183,4197,4218,4221,4224,4233,4250,4267,4296,4310,4038,4155,4162,4189,4256,4261"
list_mela=str_mela.split(",")[:] 
print(len(list_mela))
df_response_mela_raw=df_result[df_result['Type']=='Melanoma']
####################overall
print('melanoma overall 数量：',df_response_mela_raw.shape)
df_response_mela_raw['Neoantigen'] = pd.to_numeric(df_response_mela_raw['Neoantigen'], errors='coerce')  # 将非数值项转换为NaN
# 计算两两组别之间的差异
categories = df_response_mela_raw['Response'].unique()
correlations = {}
for i in range(len(categories)):
    for j in range(i + 1, len(categories)):
        cat_i = df_response_mela_raw[df_response_mela_raw['Response'] == categories[i]]['Neoantigen']
        cat_j = df_response_mela_raw[df_response_mela_raw['Response'] == categories[j]]['Neoantigen']
        u_statistic, p_value = stats.ttest_ind(cat_i, cat_j, alternative='two-sided',equal_var=False)
        correlations[(categories[i], categories[j])] = p_value
###计算响应组和不响应组之间的差异
response_list=df_response_mela_raw[df_response_mela_raw['Response']!='NR']['Neoantigen'].to_list()
noresponse_list=df_response_mela_raw[df_response_mela_raw['Response']=='NR']['Neoantigen'].to_list()
u_statistic, p_value = stats.ttest_ind(response_list, noresponse_list, alternative='two-sided',equal_var=False)
print('melenoma overall 响应Vs不响应')
print(u_statistic, p_value)
p_value=round(p_value,4)
# 绘图
# 设置绘图风格为 "ticks"（无网格线），并指定背景颜色和其他样式选项
sns.set(style="ticks", rc={"axes.facecolor": "white", "figure.facecolor": "white"})
# 定义自定义颜色
custom_palette = {'NR': '#FF9999', 'PR': '#66B2FF', 'CR': '#99FF99'}
# 创建图形和坐标轴
plt.figure(figsize=(3, 5))
# 绘制箱线图
#ax=sns.boxplot(x='Response', y='Neoantigen', data=df_response_mela, palette=custom_palette)
ax=sns.boxplot(x='Response', y='Neoantigen', data=df_response_mela_raw, boxprops = {'color':'black','facecolor':'white'},linewidth=1,medianprops = {'linestyle':'--','color':'red'})
# 添加散点（使用 stripplot 或 swarmplot）
# 使用 stripplot 可能会有重叠的点，而 swarmplot 会调整点的位置以避免重叠
#sns.stripplot(x='Response', y='Neoantigen', data=df_response_mela, palette=custom_palette, jitter=True, size=4)
sns.stripplot(x='Response', y='Neoantigen', data=df_response_mela_raw, color='black', jitter=True, size=2)
# 添加相关性注释
max_y = df_response_mela_raw['Neoantigen'].max() * 1.05  # 略高于最大值的位置放置注释
ax.annotate('Ttest(Response,non-Response)='+str(p_value), xy=(0, max_y),
                    ha='center', va='bottom', fontsize=10, color='red')
for (cat_i, cat_j), corr in correlations.items():
    if not np.isnan(corr):
        max_y-=200
        x1 = df_response_mela_raw['Response'].unique().tolist().index(cat_i)
        x2 = df_response_mela_raw['Response'].unique().tolist().index(cat_j)
        ax.annotate(f'Ttest({cat_i},{cat_j})={corr:.2f}', xy=((x1 + x2) / 2, max_y),
                    ha='center', va='bottom', fontsize=10, color='red')
# 添加标题和标签
plt.xlabel('Response')
plt.ylabel('The number of neoantigens')
plt.title('Melanoma TIL ACT :overall')
# 显示图形
outplot=args.output+"_melanoma_overall.png"
plt.savefig(outplot)
plt.show()

####################anti PD-1 experienced
df_response_mela=df_response_mela_raw[df_response_mela_raw['Sample'].isin(list_mela)]
print('melanoma anti PD-1 experienced 数量：',df_response_mela.shape)
df_response_mela['Neoantigen'] = pd.to_numeric(df_response_mela['Neoantigen'], errors='coerce')  # 将非数值项转换为NaN
# 计算两两组别之间的差异
categories = df_response_mela['Response'].unique()
correlations = {}
for i in range(len(categories)):
    for j in range(i + 1, len(categories)):
        cat_i = df_response_mela[df_response_mela['Response'] == categories[i]]['Neoantigen']
        cat_j = df_response_mela[df_response_mela['Response'] == categories[j]]['Neoantigen']
        u_statistic, p_value = stats.ttest_ind(cat_i, cat_j, alternative='two-sided',equal_var=False)
        correlations[(categories[i], categories[j])] = p_value
###计算响应组和不响应组之间的差异
response_list=df_response_mela[df_response_mela['Response']!='NR']['Neoantigen'].to_list()
noresponse_list=df_response_mela[df_response_mela['Response']=='NR']['Neoantigen'].to_list()
u_statistic, p_value = stats.ttest_ind(response_list, noresponse_list, alternative='two-sided',equal_var=False)
print('melenoma anti PD-1 experienced 响应Vs不响应')
print(u_statistic, p_value)
p_value=round(p_value,4)
# 绘图
# 设置绘图风格为 "ticks"（无网格线），并指定背景颜色和其他样式选项
sns.set(style="ticks", rc={"axes.facecolor": "white", "figure.facecolor": "white"})
# 定义自定义颜色
custom_palette = {'NR': '#FF9999', 'PR': '#66B2FF', 'CR': '#99FF99'}
# 创建图形和坐标轴
plt.figure(figsize=(3, 5))
# 绘制箱线图
#ax=sns.boxplot(x='Response', y='Neoantigen', data=df_response_mela, palette=custom_palette)
ax=sns.boxplot(x='Response', y='Neoantigen', data=df_response_mela, boxprops = {'color':'black','facecolor':'white'},linewidth=1,medianprops = {'linestyle':'--','color':'red'})
# 添加散点（使用 stripplot 或 swarmplot）
# 使用 stripplot 可能会有重叠的点，而 swarmplot 会调整点的位置以避免重叠
#sns.stripplot(x='Response', y='Neoantigen', data=df_response_mela, palette=custom_palette, jitter=True, size=4)
sns.stripplot(x='Response', y='Neoantigen', data=df_response_mela, color='black', jitter=True, size=2)
# 添加相关性注释
max_y = df_response_mela['Neoantigen'].max() * 1.05  # 略高于最大值的位置放置注释
ax.annotate('Ttest(Response,non-Response)='+str(p_value), xy=(0, max_y),
                    ha='center', va='bottom', fontsize=10, color='red')
for (cat_i, cat_j), corr in correlations.items():
    if not np.isnan(corr):
        max_y-=200
        x1 = df_response_mela['Response'].unique().tolist().index(cat_i)
        x2 = df_response_mela['Response'].unique().tolist().index(cat_j)
        ax.annotate(f'Ttest({cat_i},{cat_j})={corr:.2f}', xy=((x1 + x2) / 2, max_y),
                    ha='center', va='bottom', fontsize=10, color='red')
# 添加标题和标签
plt.xlabel('Response')
plt.ylabel('The number of neoantigens')
plt.title('Melanoma TIL ACT :anti PD-1 experienced')
# 显示图形
outplot=args.output+"_melanoma_anti_PD1_experienced.png"
plt.savefig(outplot)
plt.show()

#################anti PD-1 naïve
df_response_mela=df_response_mela_raw[~df_response_mela_raw['Sample'].isin(list_mela)]
print('melanoma anti PD-1 naïve 数量：',df_response_mela.shape)
df_response_mela['Neoantigen'] = pd.to_numeric(df_response_mela['Neoantigen'], errors='coerce')  # 将非数值项转换为NaN
# 计算两两组别之间的差异
categories = df_response_mela['Response'].unique()
correlations = {}
for i in range(len(categories)):
    for j in range(i + 1, len(categories)):
        cat_i = df_response_mela[df_response_mela['Response'] == categories[i]]['Neoantigen']
        cat_j = df_response_mela[df_response_mela['Response'] == categories[j]]['Neoantigen']
        u_statistic, p_value = stats.ttest_ind(cat_i, cat_j, alternative='two-sided',equal_var=False)
        correlations[(categories[i], categories[j])] = p_value
###计算响应组和不响应组之间的差异
response_list=df_response_mela[df_response_mela['Response']!='NR']['Neoantigen'].to_list()
noresponse_list=df_response_mela[df_response_mela['Response']=='NR']['Neoantigen'].to_list()
u_statistic, p_value = stats.ttest_ind(response_list, noresponse_list, alternative='two-sided',equal_var=False)
print('melenoma anti PD-1 naïve 响应Vs不响应')
print(u_statistic, p_value)
p_value=round(p_value,4)
# 绘图
# 设置绘图风格为 "ticks"（无网格线），并指定背景颜色和其他样式选项
sns.set(style="ticks", rc={"axes.facecolor": "white", "figure.facecolor": "white"})
# 定义自定义颜色
custom_palette = {'NR': '#FF9999', 'PR': '#66B2FF', 'CR': '#99FF99'}
# 创建图形和坐标轴
plt.figure(figsize=(3, 5))
# 绘制箱线图
#ax=sns.boxplot(x='Response', y='Neoantigen', data=df_response_mela, palette=custom_palette)
ax=sns.boxplot(x='Response', y='Neoantigen', data=df_response_mela, boxprops = {'color':'black','facecolor':'white'},linewidth=1,medianprops = {'linestyle':'--','color':'red'})
# 添加散点（使用 stripplot 或 swarmplot）
# 使用 stripplot 可能会有重叠的点，而 swarmplot 会调整点的位置以避免重叠
#sns.stripplot(x='Response', y='Neoantigen', data=df_response_mela, palette=custom_palette, jitter=True, size=4)
sns.stripplot(x='Response', y='Neoantigen', data=df_response_mela, color='black', jitter=True, size=2)
# 添加相关性注释
max_y = df_response_mela['Neoantigen'].max() * 1.05  # 略高于最大值的位置放置注释
ax.annotate('Ttest(Response,non-Response)='+str(p_value), xy=(0, max_y),
                    ha='center', va='bottom', fontsize=10, color='red')
for (cat_i, cat_j), corr in correlations.items():
    if not np.isnan(corr):
        max_y-=200
        x1 = df_response_mela['Response'].unique().tolist().index(cat_i)
        x2 = df_response_mela['Response'].unique().tolist().index(cat_j)
        ax.annotate(f'Ttest({cat_i},{cat_j})={corr:.2f}', xy=((x1 + x2) / 2, max_y),
                    ha='center', va='bottom', fontsize=10, color='red')
# 添加标题和标签
plt.xlabel('Response')
plt.ylabel('The number of neoantigens')
plt.title('Melanoma TIL ACT :anti PD-1 naïve')
# 显示图形
outplot=args.output+"_melanoma_anti_PD1_naive.png"
plt.savefig(outplot)
plt.show()

# %% [markdown]
# 绘制肺癌的新抗原和响应的关系图

# %%
df_response_lung=df_result[df_result['Type']=='Lung']
#df_response_lung=df_response_lung[df_response_lung['Response']!='SD']
df_response_lung['Neoantigen'] = pd.to_numeric(df_response_lung['Neoantigen'], errors='coerce')  # 将非数值项转换为NaN

# %%
df_response_lung['Response'].value_counts()

# %%
# 计算两两组别之间的差异
categories = df_response_lung['Response'].unique()
correlations = {}

for i in range(len(categories)):
    for j in range(i + 1, len(categories)):
        cat_i = df_response_lung[df_response_lung['Response'] == categories[i]]['Neoantigen']
        cat_j = df_response_lung[df_response_lung['Response'] == categories[j]]['Neoantigen']
        u_statistic, p_value =stats.ttest_ind(cat_i, cat_j, alternative='two-sided',equal_var=False)
        correlations[(categories[i], categories[j])] = p_value

# %%
###计算响应组和不响应组之间的差异
response_list=df_response_lung[df_response_lung['Response']!='PD']['Neoantigen'].to_list()
noresponse_list=df_response_lung[df_response_lung['Response']=='PD']['Neoantigen'].to_list()
u_statistic, p_value = stats.ttest_ind(response_list, noresponse_list, alternative='two-sided',equal_var=False)
print('Lung 响应Vs不响应')
print(u_statistic, p_value)
p_value=round(p_value,4)

# %%
# 设置绘图风格为 "ticks"（无网格线），并指定背景颜色和其他样式选项
sns.set(style="ticks", rc={"axes.facecolor": "white", "figure.facecolor": "white"})

# 定义自定义颜色
custom_palette = {'PD': '#FF9999', 'PR': '#66B2FF', 'CR': '#99FF99','SD':'blue'}

# 创建图形和坐标轴
plt.figure(figsize=(4, 5))

# 绘制箱线图
#ax=sns.boxplot(x='Response', y='Neoantigen', data=df_response_lung, palette=custom_palette)
ax=sns.boxplot(x='Response', y='Neoantigen', data=df_response_lung, boxprops = {'color':'black','facecolor':'white'},linewidth=1,medianprops = {'linestyle':'--','color':'red'})


# 添加散点（使用 stripplot 或 swarmplot）
# 使用 stripplot 可能会有重叠的点，而 swarmplot 会调整点的位置以避免重叠
#sns.stripplot(x='Response', y='Neoantigen', data=df_response_lung, palette=custom_palette, jitter=True, size=4)
sns.stripplot(x='Response', y='Neoantigen', data=df_response_lung, color='black', jitter=True, size=2)


# 添加相关性注释
max_y = df_response_lung['Neoantigen'].max() * 1.05  # 略高于最大值的位置放置注释
ax.annotate('Ttest(Response,non-Response)='+str(p_value), xy=(0, max_y),
                    ha='center', va='bottom', fontsize=10, color='red')
for (cat_i, cat_j), corr in correlations.items():
    if not np.isnan(corr):
        max_y-=100
        x1 = df_response_lung['Response'].unique().tolist().index(cat_i)
        x2 = df_response_lung['Response'].unique().tolist().index(cat_j)
        ax.annotate(f'Ttest({cat_i},{cat_j})={corr:.2f}', xy=((x1 + x2) / 2, max_y),
                    ha='center', va='bottom', fontsize=10, color='red')


# 添加标题和标签
plt.xlabel('Response')
plt.ylabel('The number of neoantigens')

# 显示图形
outplot=args.output+"_lung.png"
plt.savefig(outplot)
plt.show()


