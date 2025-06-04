###################################################################
# Created Time: Thu 25 Jan 2024 01:14:02 PM CST
#=============================================================
#!/usr/bin/env python
#-*- coding:utf8 -*-
'''
按照分组进行生存分析，输入的分组数据格式为：
This df must included: 
Harmonized_OS_Days(生存期),
Harmonized_OS_Event(生存事件),
group(分组): 分组必须写成Groupxx,
disease(癌种)
'''

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
import math
import argparse

parser = argparse.ArgumentParser(description='参数')
parser.add_argument('--input_os', '-i', help='输入样本的分组生存期')
parser.add_argument('--output', '-o', help='输出绘图结果文件')
args=parser.parse_args()


def get_exp_ci(ci_95):
	if float(ci_95)>=700:
		return round(math.exp(700),2)
	return round(math.exp(float(ci_95)),2)

def draw_lifeline_plot(df,col_name,output):
	'''This df included three columns:OS.time,OS,col_name'''
	mpl.use("agg")
	kmf = KaplanMeierFitter()
	groups = df[col_name].unique()
	#print('All',groups)
	print('All group number: ',df[col_name].value_counts())
	if len(groups)==1:
		return 0
	#print(groups,len(groups),col_name,output)
	text_num='Group Number:'
	for group in sorted(groups):
		group_data = df[df[col_name] == group]
		text_num=text_num+"\n"+group+": "+str(group_data.shape[0])
		#print('this',group_data['Harmonized_OS_Days'])
		#print('This',group,group_data.shape)
		kmf.fit(durations=group_data['Harmonized_OS_Days'], event_observed=group_data['Harmonized_OS_Event'], label=group)
		kmf.plot_survival_function()
		results = logrank_test(df['Harmonized_OS_Days'][df[col_name] == groups[0]], df['Harmonized_OS_Days'][df[col_name] == groups[1]],event_observed_A=df['Harmonized_OS_Event'][df[col_name] == groups[0]], event_observed_B=df['Harmonized_OS_Event'][df[col_name] == groups[1]])
	#plt.text(2,1,"p-value:{}".format(round(results.p_value,5)))
	#print('This',summary)
	cph=CoxPHFitter()
	cph.fit(df, 'Harmonized_OS_Days', event_col='Harmonized_OS_Event', formula=col_name)
	summary=cph.summary
	summary['HR'] = summary['exp(coef)']
	#print('This',summary)
	plt.text(10,0.18,text_num)
	plt.text(10,0.13,"p-value: {}".format(round(summary['p'].values[0],4)))
	plt.text(10,0.08,"HR: {}".format(round(summary['HR'].values[0],2)))
	plt.text(10,0.03,"95% CI: {0}-{1}".format(get_exp_ci(summary['coef lower 95%'].values[0]),get_exp_ci(summary['coef upper 95%'].values[0])))

	#if summary.p.values[0]<0.05:
		#print("###########")
		#print(col_name,output)
		#print("p-value:{}".format(round(summary.p.values[0],5)))
	print("p-value:{}".format(round(summary.p.values[0],5)))
	plt.xlabel('Time (Day)')
	plt.ylabel('OS rate')
	#plt.title('Kaplan-Meier Survival Curve')
	plt.legend(title="Neoantigen Load", loc='upper right')
	plt.title(output)
	#plt.savefig('{}.png'.format(output),dpi=300)
	#plt.savefig('{}.pdf'.format(output))
	plt.savefig(output,dpi=300)
	plt.close()

##每个癌肿分类，及其对应的数量
def check_tumor_type(df,col="disease",threshold=3):
	#print(df.head())
	counts = df[col].value_counts()
	return list(counts[counts>threshold].index)


def main():
	df=pd.read_csv(args.input_os,sep="\t")
	#print(df['group'])
	outplot=args.output
	tumor_list=check_tumor_type(df)
	print("the aim cancer: ",tumor_list)
	for k in tumor_list:
		df1=df[df["disease"]==k]
		#print(df1.columns)
		##k:癌肿，j：类别
		for j in [i for i in df1.columns if i.startswith("group")]:
			try:
				#draw_lifeline_plot(df1,j,k+"-"+j)
				#print('This',j,df1)
				draw_lifeline_plot(df1,j,outplot)
			except Exception as e:
				print("###########\t",end="")
				print(e)
if __name__=="__main__":
	main()







