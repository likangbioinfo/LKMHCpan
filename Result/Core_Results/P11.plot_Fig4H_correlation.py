# %%
import sys
import pandas as pd
import numpy as np
import re
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# %%
def creat_df(filename):
    df=pd.read_csv(filename,header=1,sep="\s+")
    file_index=df.index.tolist()
    index_list=[x[1] for x in file_index]
    index_list
    df.index =index_list
    return df

# %%
def stat_correlation(df1,df2,software):
    # 初始化一个空的列表来存储相关性值
    correlations = []
    # 计算每一行的相关性
    for i in range(len(df1)):
        row1 = df1.iloc[i].values 
        row2 = df2.iloc[i].values
        corr = np.corrcoef(row1, row2)[0, 1]
        correlations.append(round(corr,2))

    # 将相关性值转换为 DataFrame
    indexname=software
    colname='Correlation_original_'+software
    correlation_df = pd.DataFrame({indexname: df2.index, colname: correlations})

    # 输出结果
    #print(correlation_df)
    return correlation_df

# %%
def plot_fig(merged_df,title):
    # 绘制热图
    merged_df['index'] =[ str(x+1) for x in range(len(merged_df))]
    plt.figure(figsize=(6, 2))
    plt.plot(merged_df['index'], merged_df['Correlation_original_LKMHCpan2'], marker='o', linestyle='-', color='#E97123', label='LKMHCpan2')
    plt.plot(merged_df['index'], merged_df['Correlation_original_NetMHCpan4.1'], marker='o', linestyle='-', color='#E1B224', label='NetMHCpan4.1')
    #plt.xticks(merged_df['index'], merged_df['Original'])
    # 添加标题
    title=title+": Amino acid correlation for each position"
    plt.title(title)
    plt.ylabel('Correlation')
    plt.xlabel('Position')
    plt.legend(loc='lower left', fontsize=10)
    outplot="Fig_"+title+"_correlation.svg"
    plt.savefig(outplot)
    # 显示图形
    plt.show()

# %% [markdown]
# #####HLA-B*27:05

# %%
HLA2705_ori='/01.gibbscluster_HLA/HLA-B2705_C67S_3135460/matrices/gibbs.1of1.mat'
df_ori=creat_df(HLA2705_ori)

# %%
HLA2705_lk='/Table_HLAMutant_seq_data_HLA-B2705_C67S_changeAA.positive_43858/matrices/gibbs.1of1.mat'
df_lk=creat_df(HLA2705_lk)

# %%
HLA2705_net='/02.gibbscluster_predict/Table_HLAMutant_seq_data_HLA-B2705_C67S_netmhcpna4.1_83338/matrices/gibbs.1of1.mat'
df_net=creat_df(HLA2705_net)

# %%
corr_ori_lk=stat_correlation(df_ori,df_lk,'LKMHCpan2')
corr_ori_net=stat_correlation(df_ori,df_net,'NetMHCpan4.1')
joined_df= pd.concat([corr_ori_lk, corr_ori_net], axis=1)

ori_index=pd.DataFrame(df_ori.index.tolist())
ori_index.columns=['Original']
merged_df= pd.concat([ori_index, joined_df], axis=1)
merged_df.to_csv('HLA-B2705_C67S_correlation_results.csv', index=False,sep="\t")
#merged_df

# %%
plot_fig(merged_df,"HLA-B*27:05_C67S")

# %% [markdown]
# #####HLA-B*27:09

# %%
HLA2709_ori='/01.gibbscluster_HLA/HLA-B2709_C67S_3136287/matrices/gibbs.1of1.mat'
df_ori=creat_df(HLA2709_ori)

HLA2709_lk='Table_HLAMutant_seq_data_HLA-B2709_C67S_changeAA.positive_44307/matrices/gibbs.1of1.mat'
df_lk=creat_df(HLA2709_lk)

HLA2709_net='/02.gibbscluster_predict/Table_HLAMutant_seq_data_HLA-B2709_C67S_netmhcpna4.1_83361/matrices/gibbs.1of1.mat'
df_net=creat_df(HLA2709_net)

corr_ori_lk=stat_correlation(df_ori,df_lk,'LKMHCpan2')
corr_ori_net=stat_correlation(df_ori,df_net,'NetMHCpan4.1')
joined_df= pd.concat([corr_ori_lk, corr_ori_net], axis=1)

ori_index=pd.DataFrame(df_ori.index.tolist())
ori_index.columns=['Original']
merged_df= pd.concat([ori_index, joined_df], axis=1)
merged_df.to_csv('HLA-B2709_C67S_correlation_results.csv', index=False,sep="\t")
#merged_df

plot_fig(merged_df,"HLA-B*27:09_C67S")



