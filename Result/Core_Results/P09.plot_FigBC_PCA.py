# %%
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import pandas as pd
import os
from adjustText import adjust_text
import matplotlib.patches as mpatches
import argparse
import matplotlib.pyplot as plt


# %%
parser = argparse.ArgumentParser(description='参数')
parser.add_argument('--pwd', '-p', help='输入文件所在的文件夹')
parser.add_argument('--outplot', '-op', help='输出绘制的PCA聚类图')
parser.add_argument('--output', '-ow', help='输出PCA结果')
args=parser.parse_args()

# %%
filename_list=[]
# 要查找的根路径
pwd=args.pwd
# 要查找的目标路径
target_path = 'gibbs_ori_HLA'

target_path_novel = 'gibbs_HLA'

# 要查找的目标文件名
target_file = 'gibbs.1of1.mat'
merged_matrix=''
# 使用 os.walk 遍历根路径下的所有文件和文件夹
index=0
for root, dirs, files in os.walk(pwd):
    # 检查目标路径是否在当前路径中
    if target_path in root:
        for file in files:
            # 检查文件名是否为目标文件
            if file == target_file:
                file_path = os.path.join(root, file)
                #print(file_path)
                sample="Ori-"+file_path.split("/")[-3].split("_")[2].replace("HLA-","")
                #print(sample)
                filename_list.append(sample)
                index+=1
                df=pd.read_csv(file_path,skiprows=1, sep=' ')
                numpy_matrix = df.values
                #print(index,file_path)
                if index==1:
                    numpy_matrix1=numpy_matrix
                elif index==2:
                    numpy_matrix2=numpy_matrix
                    merged_matrix = np.stack([numpy_matrix1, numpy_matrix2])
                else:
                    merged_matrix = np.concatenate([merged_matrix, [numpy_matrix]], axis=0)

for root, dirs, files in os.walk(pwd):    
    if target_path_novel in root:
        for file in files:
            # 检查文件名是否为目标文件
            #print(file)
            if file == target_file:
                file_path = os.path.join(root, file)
                print(file_path)
                sample="Novel-"+file_path.split("/")[-3].split("_")[1].replace("HLA-","")
                #print(file_path.split("/")[-3])
                filename_list.append(sample)
                df=pd.read_csv(file_path,skiprows=1, sep=' ')
                numpy_matrix = df.values
                #print(index,file_path)
                merged_matrix = np.concatenate([merged_matrix, [numpy_matrix]], axis=0)
            #else:
            #    print("no file:",file)


# %%
merged_matrix.shape

# %%
# 将数据重塑为形状为 (n_samples, 180) 的二维数组
merged_matrix_data = merged_matrix.reshape(merged_matrix.shape[0], -1)

# 使用 PCA 降维
pca = PCA(n_components=2)
data_pca = pca.fit_transform(merged_matrix_data)

# 使用 K-means 聚类
kmeans = KMeans(n_clusters=3)
kmeans.fit(data_pca)

# 获取聚类结果
labels = kmeans.labels_


# %%
import colorsys
import random

num_colors = 100  # 颜色数量

# 定义颜色列表

def generate_random_color():
    r = random.random()
    g = random.random()
    b = random.random()
    return (r, g, b)

color_list = [generate_random_color() for _ in range(num_colors)]
random.shuffle(color_list)


# %%
lable_sample,mark_sample,used_color=[],[],[]
cal=0
dict_sample={}
dict_ID={}
for i in filename_list:
    if "Ori" in i:
        mark_sample.append("D")
    elif "Novel" in i:
        mark_sample.append("o")

    Id=i.replace("Ori-","").replace("Novel-","")
    if Id not in dict_sample:
        lable_sample.append(color_list[cal])
        dict_sample[Id]=color_list[cal] ##sample color
        dict_ID[color_list[cal]]=Id ##color sample
        used_color.append(color_list[cal])
        #print(Id,color_list[cal])
        #print(type(color_list[cal]))
        cal+=1
    else:
        lable_sample.append(dict_sample[Id])
        #print(Id,dict_sample[Id])


# %%
#保存PCA的结果
df_pca_wr = pd.DataFrame(data_pca)
df_list = pd.DataFrame(filename_list)
df_final = pd.concat([df_list, df_pca_wr], axis=1)
df_final.columns=['Sample','PCA1','PCA2']
df_final.to_csv(args.output, index=False)

# %%
fig = plt.figure(figsize=(5, 4))
dpi = 800
# 绘制散点图
sizes=230
for i, (x, y) in enumerate(data_pca):
    #print(x,y)
    plt.scatter(x, y, c=lable_sample[i], cmap='jet',s=sizes,marker=mark_sample[i],alpha=0.8)


categories = lable_sample
colors = used_color
# 创建图例
legend_elements = []
for category, color in zip(categories, colors):
    ID=dict_ID[color]
    legend_elements.append(mpatches.Patch(color=color, label=f'HLA- {ID}'))

# 设置图例标题
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('Principal Component Analysis')

# 显示图形
plt.savefig(args.outplot, dpi=dpi,bbox_inches='tight')
plt.show()


