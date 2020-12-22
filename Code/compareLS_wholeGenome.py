#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 12:07:52 2020

@author: lu
"""
import statistics as ss
import readCSVfiles as RCF
import numpy as np
import matplotlib.pyplot as plt
from math import log
from matplotlib import colors as mcolors
import seaborn as sns
from matplotlib.font_manager import FontProperties


fname="../csvs/6LS_Selection.csv"


matrix,rowHeaderList,columnHeaderList= RCF.readCSV(fname,readAsString=True)


colorSet=(sns.color_palette())
print(columnHeaderList)

sortedMatrix=[]

wholeGenomeIndex=7
wholeGenomeList=matrix[wholeGenomeIndex]
wholeGenomeList=[float(i) for i in wholeGenomeList]
wholeGenomeList=np.asarray(wholeGenomeList)
idx= np.asarray(wholeGenomeList).argsort()

wholeGenomeList=wholeGenomeList[idx]


for lst in matrix:
    lst=np.array(lst)
    sortedMatrix.append(lst[idx])
    
codonName=sortedMatrix[0]

plt.figure()
import matplotlib.colors as colors
sorted_color_names = list(colors._colors_full_map.values())



rowHeaderList=rowHeaderList[1:]

x=range(0,len(codonName))
for i in range(1,len(sortedMatrix)-1):
    print(rowHeaderList[i])
    lst=sortedMatrix[i]
    lst=[float(i) for i in lst]
    plt.plot(x,lst,c=colorSet[i],label=rowHeaderList[i])

plt.plot(x,wholeGenomeList,label="Whole Genome",c="grey")
    
fontP = FontProperties()
fontP.set_size('small')

plt.xticks(x, codonName, rotation='vertical',size=5)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=4, prop={'size': 6})
plt.savefig("selection_comparison.pdf",bbox_inches='tight')
        
    
    

#print(matrix)