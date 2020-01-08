##!/usr/bin/env python2
## -*- coding: utf-8 -*-
#"""
#Created on Fri Oct 25 11:16:41 2019
#
#@author: lu
#"""
#
from numpy.polynomial.polynomial import polyfit

import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale
import findSequenceById
import windowPhi as methods
import numpy as np
import scipy.stats as ss
import math
import random as rd
import statsmodels.api as sm
import time
rd.seed(0)
tissueSamples=25

methods.deltaEtaFile="/home/lu/Desktop/data_modified/Crei_Selection.csv"
methods.deltaMFile="/home/lu/Desktop/data_modified/Crei_Mutation.csv"
sequenceDict=findSequenceById.findSequenceByID("/home/lu/Desktop/sequences/c_elegan.fasta",idType="gene")



f=open("cEl.csv","r+")
lines=f.readlines()

#21 tissues or lifestages

tissueMatrix=[]
phiList=[]

for i in range(tissueSamples):
    tissueMatrix.append([])

for line in lines:
    randomN=rd.randrange(1,101)
    if randomN>=5:
        continue
    line=line.rstrip()
    splitList=line.split(",")
    geneName=splitList[0]
    if geneName in sequenceDict:
        sequence=sequenceDict[geneName]
        codonList=methods.loadSequence(sequence)
        if not len(codonList)==0:
            phi=methods.calPhiForGene(sequence)
        if math.isnan(phi):
            phiList.append(0)  
        else:
            phiList.append(phi)
    else:
        continue
    for i in range(1,len(splitList)):
        if len(splitList[i])==0:
            tissueMatrix[i-1].append(0)
        else:
            tissueMatrix[i-1].append(float(splitList[i]))



def corelate(list1,list2):
    return ss.pearsonr(list1,list2)[0]

for tissueList in tissueMatrix:
    print ss.pearsonr(tissueList,phiList)

larva_L1=tissueMatrix[7]
adult=tissueMatrix[12]
embryos=tissueMatrix[16]
corelationList=[]



#print corelate(larva_L1,phiList)
#print corelate(adult,phiList)
#print max(corelationList)
    

def calculateWeight(list1,list2):
    A = np.vstack([list1, np.ones(len(list1))]).T
    m, c = np.linalg.lstsq(A,list2, rcond=None)[0]
    return m

#def regress(list1,list2):
#    b=list2
#    a=list1
#    w= np.random.uniform(1, 20, size=(899,))
#    b=sm.add_constant(b)
#    mod_wls = sm.WLS(a, b, weights=w)
#    res_wls = mod_wls.fit()
#    print(res_wls.summary())
    
def getCorelation(a,b):
    return ss.pearsonr(a,b)[0]

def testWeights(list1,list2,phiList):
    maxCorelation=0
    maxI=0
    maxJ=0
    weights=[0,1,2,3,4,5,6,7,8,9,10,100]
    for i in weights:
        for j in weights:
             newList=[]
             for k in range(len(list1)):
                 newW= i* list1[k] + j* list2[k]
                 newList.append(newW)
             cor=float(ss.pearsonr(newList,phiList)[0])
             if cor > maxCorelation:
                 maxCorelation=cor
                 maxI=i
                 maxJ=j
#             print time.time()-startTime
    print "max I: %f"%(maxI)
    print "max J: %f"%(maxJ)
    print "max Corelation: %f"%(maxCorelation)
                    


def getWeights(list1,list2,phiList):
    list1=minmax_scale(list1)
    list2=minmax_scale(list2)
    phiList=minmax_scale(phiList)
    
    weights=[-1,1,2,3,4,5]
    minDiff=99999
    maxI=0
    maxJ=0
    for i in weights:
        for j in weights:
             newList=[]
             for k in range(len(list1)):
                 diff= abs(i* list1[k] + j* list2[k] - phiList[k])
                 newList.append(diff)
             sumDiff=sum(newList)
             if sumDiff< minDiff:
                 minDiff=sumDiff
                 maxI=i
                 maxJ=j
#    print maxI
#    print maxJ
#    print minDiff



a= np.array(phiList)
#print phiList
b= np.array(adult)


#
#import numpy as np
#from sklearn import preprocessing
##print "=="
##print ss.pearsonr(embryos,adult)
##regress(embryos,adult)
#getWeights(embryos,adult,phiList)
#print getCorelation(embryos,phiList)
#print getCorelation(adult,phiList)
#
#getWeights(larva_L1,adult,phiList)
#print getCorelation(larva_L1,phiList)
#print getCorelation(adult,phiList)



def logify(lst):
    logList=[]
    for i in lst:
        if i==0:
            i=0.00001
        logList.append(np.log(i))
    return logList
    

a= np.array(phiList)
#print phiList
b= np.array(adult)

def regress(list1,list2,weightList=[]):
    if len(weightList)==0:
        model = sm.OLS(list2, list1).fit()
    else:
        model = sm.WLS(list2, list1, weights=weightList).fit()
    predictions = model.predict(list1)
    print model.params
    print model.summary()
    return predictions
#Below is a test use case for regress function, Y is the dependent factor and X are the features, the first column in X is intercept
#Y = [1,3,4,5,2,3]
#X = [[0,1,2],[0,2,3],[0,3,3],[0,5,5],[0,3,1],[0,4,4]]
#print regress(X,Y)

indexList=[3,4,5,6,8,9,11,15,18,19,20,21,22,23,24,25,26]


tempMatrix=[]
for i in range(len(tissueMatrix)):
    if i in indexList:
        tempMatrix.append(tissueMatrix[i])
        
tissueMatrix=tempMatrix


tissueMatrix_norm=[]
for tissueList in tissueMatrix:
    tissueMatrix_norm.append((minmax_scale(tissueList)))



newIndexList=[]
for i in indexList:
    newIndexList.append(i-1)

tissueMatrix_norm_trans=np.transpose(tissueMatrix_norm)
print "*****"
print tissueMatrix_norm_trans.shape 
phiList=minmax_scale(phiList)

X=sm.add_constant(tissueMatrix_norm_trans)

prediction= regress(X,phiList)

b,m=polyfit(phiList,prediction,1)
plt.plot(phiList, b + m * phiList, '-')
plt.show()

for tissueList in tissueMatrix:
    print corelate(phiList,tissueList)
print "=="
print corelate(prediction,phiList)

import numpy as np
import matplotlib.pyplot as plt

# Initiate some data, giving some randomness using random.random().
x = tissueMatrix_norm[0]
y = phiList




