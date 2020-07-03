#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 14:13:07 2020

@author: lu
"""

import numpy as np; 
import readCSVfiles as RCF
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot
import scipy.stats as ss

def testCorelation(x,y,corelationFunction):
    if "pearson" in corelationFunction:
#        print ("p value %f"%stats.pearsonr(x, y)[1])
        return  ss.pearsonr(x, y)[0]
    elif "spearman" in corelationFunction:
        return ss.spearmanr(x,y,nan_policy="omit")[0]
    elif "kendall" in corelationFunction:
        return ss.kendalltau(x,y,nan_policy="omit")[0]


#This creates a figure consisits of n*(n-1)/2 subplots, each is a 1v1 corelation
def plotHeatMap(matrix,header,savefile=False):
    import seaborn as sns; sns.set(color_codes=True)
    plt.figure(figsize=(15, 15))
    matrix=np.transpose(matrix)
#    print(rowHeader)
    g = sns.heatmap(matrix, cmap='hot',vmin=0, vmax=100)
    g.set_xlabel("Life Stages")
    g.set_ylabel("Gene")
    if savefile:
        plt.savefig("heatmap.pdf")
    else:
        plt.show()
    
    
    
    

matrix,rowHeader,columnHeader=RCF.readCSV("/home/lu/AcrossTissue/csvs/cEl.csv")

#plotHeatMap(matrix,rowHeader)
def calculateCorelation1V1(matrix,header,corelationMethod="pearson"):
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i<j:
                ls1=matrix[i]
                ls2=matrix[j]
                r=testCorelation(ls1,ls2,corelationFunction=corelationMethod)
#                print (header[i])
calculateCorelation1V1(matrix,rowHeader)
print(len(rowHeader))

    