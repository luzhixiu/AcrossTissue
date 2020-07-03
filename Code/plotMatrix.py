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


#This creates a figure consisits of n*(n-1)/2 subplots, each is a 1v1 corelation
def plotHeatMap(matrix,header,savefile=False):
    import seaborn as sns; sns.set(color_codes=True)
    plt.figure(figsize=(15, 15))
    matrix=np.transpose(matrix)
    print(rowHeader)
    g = sns.heatmap(matrix, cmap='hot',vmin=0, vmax=100)
    g.set_xlabel("Life Stages")
    g.set_ylabel("Gene")
    if savefile:
        plt.savefig("heatmap.pdf")
    else:
        plt.show()
    
    
    
    

matrix,rowHeader,columnHeader=RCF.readCSV("/home/lu/AcrossTissue/csvs/cEl.csv")

plotHeatMap(matrix,rowHeader)
