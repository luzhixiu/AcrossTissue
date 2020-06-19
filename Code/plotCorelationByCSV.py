#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 11:40:32 2020

@author: lu
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
inputFile="/home/lu/AcrossTissue/ComparisonLSG/selectionParameterBetweenLS-source.csv"
headers=[]


#format: first line should be the header, each column is a list
def parseInputFile(fname):
    csvMatrix=[]
    f=open(fname)
    lines=f.readlines()
    global headers
    print lines[0]
    headers=lines[0].rstrip().split(",")
    for line in lines[1:]:
        splitList=line.rstrip().split(",")
        print splitList
        splitList = list(map(float, splitList)) 
        csvMatrix.append(splitList)
    return csvMatrix

#====================private helper functions
def estimate_coef(x, y): 
    # number of observations/points 
    n = np.size(x) 
    # mean of x and y vector 
    m_x, m_y = np.mean(x), np.mean(y) 
    # calculating cross-deviation and deviation about x 
    SS_xy = np.sum(y*x) - n*m_y*m_x 
    SS_xx = np.sum(x*x) - n*m_x*m_x 
    # calculating regression coefficients 
    b_1 = SS_xy / SS_xx 
    b_0 = m_y - b_1*m_x 
    return(b_0, b_1)  
    
def testCorelation(x,y,corelationFunction):
    if "pearson" in corelationFunction:
        print ("p value %f"%ss.pearsonr(x, y)[1])
        return  ss.pearsonr(x, y)[0]
    elif "spearman" in corelationFunction:
        return ss.spearmanr(x,y,nan_policy="omit")[0]
    elif "kendall" in corelationFunction:
        return ss.kendalltau(x,y,nan_policy="omit")[0]    
    
    

    
def plotCorelation(x,y,a,b,xLabel="list1",yLabel="list2",logScale="no",showCorelation="yes"):
    global ax
    subAxis=ax[a,b]
    if(logScale=="yes"):
        x=logify(x)
        y=logify(y)
    if "yes" in showCorelation:
        subAxis.text(0.1,0.68,"R =: %0.4f"%(testCorelation(x,y,"pearson")))
    x=np.array(x)
    y=np.array(y)
    subAxis.set_xlabel(xLabel)
    subAxis.set_ylabel(yLabel)

    b=estimate_coef(x,y)
    y_pred=b[0]+b[1]*x
    subAxis.plot(x, y_pred, color = "g",linewidth=1)
    subAxis.scatter(x,y,s=0.5)


csvMatrix=parseInputFile(inputFile)

matrix=np.transpose(np.asarray(csvMatrix))



plt.figure()
n=len(headers)

fig, ax = plt.subplots(n, n,figsize=(20, 12))
fig.tight_layout(pad=1.5)


#ax[0, 0].plot(range(10), 'r') #row=0, col=0
#ax[1, 0].plot(range(10), 'b') #row=1, col=0
#ax[0, 1].plot(range(10), 'g') #row=0, col=1
#ax[1, 1].plot(range(10), 'k') #row=1, col=1


print headers
for a in range(len(headers)):
    for b in range(len(headers)):
        print matrix[a]
        print testCorelation(matrix[a],matrix[b],"pearson")
        plotCorelation(matrix[a],matrix[b],a,b,xLabel=headers[a],yLabel=headers[b])
plt.show()
fig.savefig(inputFile+"_plot.pdf",bbox_inches='tight')




























