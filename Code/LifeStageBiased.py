#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:18:34 2020

@author: lu
"""
import statistics as ss
import readCSVfiles as RCF
import numpy as np
import matplotlib.pyplot as plt
from math import log


#A few preset options here:

# Cut off of low exp (less than median) switch


def logify(lst):
    return [log(y,10) for y in lst]

#Output format should look like this: {GeneId,LS,LS_exp,secondMax,restMean,foldDiff}
def findLSGene(expMatrix,rowHeaderList,columnHeaderList,foldDiffcutOff=2):
    outputFinal=""
    geneIdList=[]
    LSList=[]
    expList=[]
    secondMaxList=[]
    restMeanList=[]
    foldDiffList=[]
    
    headerStr="GeneID,LS,LS_EXP,SecondMax,RestMean,FoldDiff \n"
    outputFinal+=headerStr
    global expMedian
    if not cutLowExp:
        expMedian=0 #this removes the low cutoff of expression, comment this line to filter out lowely expressed genes.
    LS_Genes_Count=0
    for k in range(len(expMatrix)):
        ls=expMatrix[k]
        for i in range(len(ls)):
            exp=ls[i]
            if exp>= expMedian:
                restMean= (sum(ls)-exp)/(len(ls)-1)
                secondMax=sorted(ls,reverse=True)[1]
                if exp>restMean*foldDiffcutOff:
                    foldDiff=exp/restMean
                    geneId=rowHeaderList[k]
                    LS=columnHeaderList[i]
                    LS_Genes_Count+=1
                    geneIdList.append(geneId)
                    LSList.append(LS)
                    expList.append(exp)
                    secondMaxList.append(secondMax)
                    restMeanList.append(restMean)
                    foldDiffList.append(foldDiff)
                    outputString=geneId+","+ str(LS) + "," + str(exp) + "," + str(secondMax)+ "," + str(restMean) + "," + str(foldDiff)+ "\n" 
                    outputFinal+=outputString
    
#    expList=logify(expList)
#    foldDiffList=logify(expList)
    return outputFinal
        
def main():
    cutLowExp=False
    matrix,rowHeaderList,columnHeaderList= RCF.readCSV("/home/lu/AcrossTissue/csvs/5LS_L2L3Combined.csv")
    
    LSnames=rowHeaderList[1:]
    matrix=np.transpose(matrix)
    
    
    expSum=0.0
    expCount=0
    expList=[]
    for ls in matrix:
        for exp in ls:
            if exp==0:
                continue
            expSum+=exp
            expCount+=1
            expList.append(exp)
    expMean=expSum/expCount
    expMedian=ss.mean(expList)

    assert len(columnHeaderList),len(matrix)
    outputS=findLSGene(matrix,columnHeaderList,LSnames)
    f=open("/home/lu/AcrossTissue/csvs/LifeStageGenes_collapse.csv","w")
    f.write(outputS)
    f.close()
    