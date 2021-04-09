from pathlib import Path
import os
import pandas as pd
import scipy.stats as ss
import seaborn as sns
import matplotlib.pyplot as plt
inputDir="/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/Top_300_From_Each_Group"
inputDir="/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/Fold_Diff_4"


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
    df_sub=df_LS[["WBID","PHI","Std.Dev"]]
    nameConv={'PHI':("PHI_"+name),"Std.Dev":"StdDev_"+name}
    corTarget.append("PHI_"+name)
    df_sub=df_sub.rename(columns=nameConv)
    df=pd.merge(df,df_sub,on="WBID",how="outer")


df=pd.merge(df_empirical,df,on="WBID",how="outer")

corTarget.append("Phi_WholeGenome")
target=["adult Ce","PHI_adult","elongating embryo Ce",]
print(corTarget)
df_cor=df[corTarget]

df_cor=df_cor.sort_index(axis=1)


cor=df_cor.corr(method="spearman")

sns.heatmap(cor, annot=True,cmap="magma_r",annot_kws={"size": 8})
# plt.xticks(rotation=-25,size=5) 
# plt.xlabel('xlabel', fontsize=1)
print(cor)
