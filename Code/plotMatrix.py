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


#This creates a figure consisits of n*(n-1)/2 subplots, each is a 1v1 corelation
def plotCorelation_1v1(matrix,header):
    plt.figure()
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.show()
    
    
    
    

matrix,rowHeader,columnHeader=RCF.readCSV("/home/lu/AcrossTissue/csvs/cEl.csv")
print matrix
print rowHeader
#print columnHeader
plotCorelation_1v1(np.asarray(matrix),header)


