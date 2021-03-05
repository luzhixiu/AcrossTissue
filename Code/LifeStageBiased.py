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
inputFile=""
outputFile=""
cutLowExp=True 
cutLowPercentile=0.1
fixedCutValue=0
    
#In code global parameter   
matrix,rowHeaderList,columnHeaderList=[],[],[]


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
    
    headerStr="GeneID,LS,LS_EXP,SecondMax,RestMean,FoldDiff\n"
    outputFinal+=headerStr
    global cutLowExp,cutLowPercentile
    cutOffValue=0
    if cutLowExp:
        #this removes the low cutoff of expression, comment this line to filter out lowely expressed genes.
        flatMatrix=sorted(np.array((expMatrix).flatten()))
        cutOffValue=((flatMatrix)[int(len(flatMatrix)*cutLowPercentile)])
    if fixedCutValue>0:
       cutOffValue=fixedCutValue 
    LS_Genes_Count=0
    for k in range(len(expMatrix)):
        ls=expMatrix[k]
        for i in range(len(ls)):
            exp=ls[i]
            if exp>= cutOffValue:
                restMean= (sum(ls)-exp)/(len(ls)-1)
                secondMax=sorted(ls,reverse=True)[1]
                if exp>secondMax*foldDiffcutOff:
                    if secondMax==0:
                        foldDiff=float('inf')
                    else:
                        foldDiff=exp/secondMax
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



def readFile():
    global matrix, rowHeaderList,columnHeaderList,inputFile
    matrix,rowHeaderList,columnHeaderList= RCF.readCSV(inputFile)
    
        
def main():
    global matrix,rowHeaderList,columnHeaderList
    if len(matrix)<=1:
        readFile()
    LSnames=rowHeaderList[1:]
    matrix=np.transpose(matrix)
    assert len(columnHeaderList),len(matrix)
    outputS=findLSGene(matrix,columnHeaderList,LSnames)
    global outputFile
    f=open(outputFile,"w")
    f.write(outputS)
    f.close()
