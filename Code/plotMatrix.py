#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 14:13:07 2020

@author: lu
"""

import numpy as np;
import readCSVfiles as RCF
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot
import pandas as pd
import scipy.stats as ss

np.set_printoptions(precision=2)
matplotlib.use('Qt5Agg')

def testCorelation(x, y, corelationFunction):
    if "pearson" in corelationFunction:
        #        print ("p value %f"%stats.pearsonr(x, y)[1])
        return ss.pearsonr(x, y)[0]
    elif "spearman" in corelationFunction:
        return ss.spearmanr(x, y, nan_policy="omit")[0]
    elif "kendall" in corelationFunction:
        return ss.kendalltau(x, y, nan_policy="omit")[0]


# This creates a figure consisits of n*(n-1)/2 subplots, each is a 1v1 corelation
def plotHeatMap(matrix, header, savefile=False):
    import seaborn as sns;
    sns.set(color_codes=True)
    plt.figure(figsize=(15, 15))
    matrix = np.transpose(matrix)
    #    print(rowHeader)
    g = sns.heatmap(matrix, cmap='hot', vmin=0, vmax=100)
    g.set_xlabel("Life Stages")
    g.set_ylabel("Gene")
    if savefile:
        plt.savefig("heatmap.pdf")
    else:
        plt.show()


# plotHeatMap(matrix,rowHeader)
def calculateCorelation1V1(matrix, header, corelationMethod="pearson"):
    corelationMatirx = []
    for i in range(len(matrix)):
        corelationList = []
        for j in range(len(matrix)):
            ls1 = matrix[i]
            ls2 = matrix[j]
            r = testCorelation(ls1, ls2, corelationFunction=corelationMethod)
            corelationList.append(round(r, 4))
        corelationMatirx.append(corelationList)
    #    print (corelationMatirx)
    corelationMatirx = np.array(corelationMatirx)
    df = pd.DataFrame(corelationMatirx, columns=header)
    fig, ax = plt.subplots(dpi=200)
    # hide axes
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    tb = ax.table(cellText=df.values,colLabels=header,rowLabels=header,loc='center')
    tb.auto_set_font_size(False)
    tb.set_fontsize(2.8)
    # fig.tight_layout()
    fig.savefig("../Graph/allLSCorelation.pdf")
    plt.show()

def plotClusterMap(matrix,rowHeader):
    import seaborn as sns;
    matrix=np.transpose(matrix)
    df=pd.DataFrame(matrix, columns =rowHeader)
    sns.set(color_codes=True)
    g = sns.clustermap(df,vmin=0, vmax=100)
    plt.savefig("../Graph/ClusterMap.pdf")

matrix, rowHeader, columnHeader = RCF.readCSV("/home/lu/AcrossTissue/csvs/cel_LS.csv")
# remove all the genes that have at least one unobserved value in all these life stages
matrix=np.transpose(matrix)
fullMatrix=[]
for ls in matrix:
    if not 0 in ls:
        fullMatrix.append(ls)
matrix=np.transpose(fullMatrix)

# =========================================================================================================

# matrix,rowHeader,columnHeader=RCF.readCSV("/home/lu/AcrossTissue/csvs/cEl.csv",partialIndexList=[5,6,7,8,9])
# rowHeader = [header.replace("hermaphrodite__organism__", "") for header in rowHeader]
# rowHeader = [header.replace("_Ce", "") for header in rowHeader]
# rowHeader = [header.replace("_hermaphrodite", "") for header in rowHeader]
# print (rowHeader)
# print (np.shape(matrix))



# calculateCorelation1V1(matrix, rowHeader)

rowHeader=rowHeader[1:]
plotClusterMap(matrix,rowHeader)