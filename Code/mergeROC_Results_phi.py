from pathlib import Path
import os
import pandas as pd
import scipy.stats as ss
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def r2(x, y):
    return ss.pearsonr(x, y)[0] ** 2


# inputDir="/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/Top_300_From_Each_Group"
inputDir="/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/Fold_Diff_4"
# 

wholeGenomePhiPath='/home/lu/AcrossTissue/csvs/wholeGenomePhi_WBID.csv'
empirical_LS_Path='/home/lu/AcrossTissue/csvs/5LS_empirical_exp.csv'

rootdir = inputDir


phi_result_dir_list=[]
fasta_file_dir_list=[]
def listdirs(rootdir):
    global phi_result_dir_list
    for it in os.scandir(rootdir):
        if it.is_dir():
            listdirs(it)
        elif("gene_expression" in str(it)):
            phi_result_dir_list.append(it.path)
            
listdirs(rootdir)

for dr in phi_result_dir_list:
    parentDir=Path(dr).parents[2]
    for it in os.scandir(parentDir):
        p=it.path
        if "fasta" in p:
            fasta_file_dir_list.append(p)


corTarget=[]            
df_empirical=pd.read_csv(empirical_LS_Path) 
df_empirical=df_empirical.add_prefix("EMP_")
print(df_empirical.columns)
df_empirical=df_empirical.rename(columns={"EMP_WBID":"WBID"})


corTarget.extend(list(df_empirical.columns))
df = pd.read_csv(wholeGenomePhiPath)    


for i in range(len(fasta_file_dir_list)):
    fastaPath=fasta_file_dir_list[i]
    f=open(fastaPath,"r")
    name=(fastaPath.split("/")[-2])
    lines=f.readlines()
    wormbaseIdList=[]
    for line in lines:
        if ('>' in line):
            wormbaseId=line[1:].rstrip()
            wormbaseIdList.append(wormbaseId)
    
    phi_result_dir=phi_result_dir_list[i]
    df_LS=pd.read_csv(phi_result_dir)
    df_LS.insert(loc=0,column="WBID",value=wormbaseIdList)
    df_sub=df_LS[["WBID","Mean","Std.Dev"]]
    nameConv={'Mean':("PHI_"+name),"Std.Dev":"StdDev_"+name}
    corTarget.append("PHI_"+name)
    df_sub=df_sub.rename(columns=nameConv)
    df=pd.merge(df,df_sub,on="WBID",how="outer")


df=pd.merge(df_empirical,df,on="WBID",how="outer")

corTarget.append("Phi_WholeGenome")

target=[]
colNames=list(df.columns)
for colName in colNames:
    if "PHI" in colName or "EMP" in colName:
        target.append(colName)
print(target)



df_log=pd.DataFrame()
numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
for c in [c for c in df.columns if df[c].dtype in numerics]:
    df_log[c] = np.log10(df[c]+1)


df_cor=df_log[target]
df_cor=df_cor.sort_index(axis=1)

cor=df_cor.corr(method="pearson")


cor_x=[]
for x in cor.columns:
    if "EMP" in x:
        cor_x.append(x)

cor=cor[cor_x]
cor=cor.filter(like="PHI",axis=0)
df_log=df_log.sort_index(axis=1,ascending=False)


plt.figure(dpi=300)
sns.heatmap(cor, annot=True,cmap="magma_r",annot_kws={"size": 8})
plt.xticks(rotation=-25) 
# plt.xlabel('xlabel', fontsize=1)

g=sns.pairplot(vars=target, data=df_log,kind="reg",corner=True,plot_kws={'line_kws':{'color':'black'}})
g.fig.suptitle("log_exp")