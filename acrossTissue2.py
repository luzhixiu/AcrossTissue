#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 13:28:36 2019

@author: lu
"""

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
import numpy as np
import scipy.stats as ss
import math
import random as rd
import statsmodels.api as sm
import time
import sys

rd.seed(0)
tissueSamples=25

def write_list_to_file(guest_list, filename):
    """Write the list to csv file."""

    with open(filename, "w") as outfile:
        for entries in guest_list:
            outfile.write(entries)
            outfile.write("\n")


def getSequencialTags(inputFile,idType="locus_tag"):
    print ("Selected id Type: %s"%(idType))
    from Bio import SeqIO
    records=SeqIO.parse(inputFile, "fasta")
    cnt=0
    mySum=0
    tagList=[]
    for record in records:
        mySum+=1
        header=str(record.description)
        startTargetIndex=header.find(str(idType))
        
        if startTargetIndex<0:
#            print "couldn't find the target idType"
            cnt+=1
            continue
        startIndex=startTargetIndex+len(idType)+1
        idName=""
        charIndex=startIndex
        while not header[charIndex]=="]":
            idName+=header[charIndex]
            charIndex+=1
        tagList.append(idName)
    return tagList

def loadSequence(sequence):
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    codonList=[]
    i=0
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            codonList.append(codon)
        i+=3
    actualCodonList=[]
    started=False
    for codon in codonList:
        if codon in stopCodonList:
            break
        if started:
            actualCodonList.append(codon)
        if codon==startCodon:
            started=True
    codonList=actualCodonList
   # print "codon readed successful, the number of codon in this sequence is %d"%(len(codonList))
    return codonList


def logify(lst):
    logList=[]
    for i in lst:
        if i==0:
            i=0.00001
        logList.append(np.log(i))
    return logList
 
    
def regress(list1,list2,weightList=[]):
    if len(weightList)==0:
        print "using OLS"
        model = sm.OLS(list2, list1).fit()
    else:
        print "using WLS"
        model = sm.WLS(list2, list1, weights=weightList).fit()
    predictions = model.predict(list1)
    print model.params
    print model.summary()
    return predictions

geneTagList=[]
phiList=[]
stdEList=[]

f=open("gene_expression.csv","r+")
lines=f.readlines()
header=lines[0]
lines=lines[1:]
geneEstDict=dict()


for line in lines:
#    print line
    
    splitList=line.split(",")
    geneTagList.append(splitList[0]) #logPhi is at 2 and log stdE is at 4
    geneTag=splitList[0]
    phi=float(splitList[1])
    stdE=float(splitList[4])
    phiList.append(phi)
    stdEList.append(stdE)
    geneEstDict[geneTag]=[phi,stdE]

#print geneEstDict
print geneEstDict["homt-1"]

f=open("cEl.csv","r+")
lines=f.readlines()

#21 tissues or lifestages

allMeasurement=[]


tissueMatrix=[]
phiList=[]
indexList=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

print len(indexList)
stdEList=[]

for i in range(tissueSamples):
    tissueMatrix.append([])
cnt=0
for line in lines:
    line=line.rstrip()
    splitList=line.split(",")
    geneName=splitList[0]
    if geneName in geneEstDict:
        phi=float(geneEstDict[geneName][0])
        stdE=geneEstDict[geneName][1]
        phiList.append(phi)
        stdEList.append(stdE)
    else:
        continue
    for i in range(1,len(splitList)):
        marked=False
        if len(splitList[i])==0:
            tissueMatrix[i-1].append(0)
            marked=True
            
        else:
            tissueMatrix[i-1].append(float(splitList[i]))
            allMeasurement.append(float(splitList[i]))
        if marked:
            cnt+=1
print marked
print min(allMeasurement)
print max(allMeasurement)
print "-----------------"

def corelate(list1,list2):
    return ss.pearsonr(list1,list2)[0]

sampleLables=["hermaphrodite, NSM, L1 larva Ce","hermaphrodite, gonad, adult Ce", "hermaphrodite, motor neuron, L2 larva Ce", "hermaphrodite, neuron, L1 larva Ce" ,"hermaphrodite, organism, 3-fold embryo Ce","hermaphrodite, organism, 4-cell embryo Ce","hermaphrodite, organism, L1 larva Ce","hermaphrodite, organism, L2 larva Ce","hermaphrodite, organism, L2d-dauer molt Ce","hermaphrodite, organism, L3 larva Ce","hermaphrodite, organism, L4 larva Ce","hermaphrodite, organism, adult Ce","hermaphrodite, organism, dauer larva Ce","hermaphrodite, organism, elongating embryo Ce","hermaphrodite, organism, enclosing embryo Ce","hermaphrodite, organism, fully-elongated embryo Ce","hermaphrodite, organism, gastrulating embryo Ce","hermaphrodite, organism, late cleavage stage embryo Ce","hermaphrodite, organism, newly molted young adult hermaphrodite Ce","hermaphrodite, organism, post dauer stage Ce","hermaphrodite, organism, proliferating embryo Ce","hermaphrodite, pharyngeal muscle cell, fully-elongated embryo Ce","hermaphrodite, somatic cell, embryo Ce","male, organism, L4 larva Ce","male, organism, embryo Ce"
]


print len(sampleLables)
print len(tissueMatrix)

print tissueMatrix[0][0]




for i in range(len(tissueMatrix)):
    if i in indexList:
        tissueList=tissueMatrix[i]
        label=sampleLables[i]
        print label
        print ss.pearsonr(tissueList,phiList)
        plt.figure()
        plt.title(label)
        plt.hist((tissueList),bins=10,log=True)
        plt.show()

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


#Below is a test use case for regress function, Y is the dependent factor and X are the features, the first column in X is intercept
#Y = [1,3,4,5,2,3]
#X = [[0,1,2],[0,2,3],[0,3,3],[0,5,5],[0,3,1],[0,4,4]]
#print regress(X,Y)


print len(tissueMatrix)

tempMatrix=[]
for i in range(len(tissueMatrix)):
    if i in indexList:
        tempMatrix.append(tissueMatrix[i])

print len(tempMatrix)
print "--"
tissueMatrix=tempMatrix


tissueMatrix_norm=[]
for tissueList in tissueMatrix:
    tissueMatrix_norm.append(minmax_scale(tissueList))

for i in range(len(tissueMatrix_norm)):
    if i in indexList:
        tissueList=tissueMatrix[i]
        label=sampleLables[i]
        print label
        print ss.pearsonr(tissueList,phiList)
        plt.figure()
        plt.title(label)
        plt.hist((tissueList),bins=10,log=True)
        plt.show()

newIndexList=[]
for i in indexList:
    newIndexList.append(i-1)

tissueMatrix_norm_trans=np.transpose(tissueMatrix_norm)
print "*****"
print tissueMatrix_norm_trans.shape 
#phiList=minmax_scale(phiList)

X=sm.add_constant(tissueMatrix_norm_trans)

prediction= regress(X,phiList)


#b,m=polyfit(phiList,prediction,1)
#plt.plot(phiList, b + m * phiList, '-')
#plt.show()

phiList=minmax_scale(phiList)


for tissueList in tissueMatrix:
    print corelate(phiList,tissueList)
    
print "=="
print corelate(prediction,phiList)

import numpy as np
import matplotlib.pyplot as plt

# Initiate some data, giving some randomness using random.random().
x = tissueMatrix_norm[0]
y = phiList

prediction= regress(X,phiList,weightList=stdEList)


#b,m=polyfit(phiList,prediction,1)
#plt.plot(phiList, b + m * phiList, '-')
#plt.show()

print "=="
print corelate(prediction,phiList)

import numpy as np
import matplotlib.pyplot as plt

# Initiate some data, giving some randomness using random.random().
x = tissueMatrix_norm[0]
y = phiList

import validator
validator.validate(phiList,prediction,logScale="no",xLabel="Phi",yLabel="Prediction")

from sklearn.decomposition import PCA
def makePCA(X,n=5):
    pca = PCA(n_components=n)
    pca.fit(X)
    PCA(n_components=n)
    pcaRatio=pca.explained_variance_ratio_
    plt.figure()
    labels=["PC1","PC2","PC3","PC4","PC5"]
    plt.bar(x=range(1,6),height=pcaRatio,tick_label=labels)
    plt.ylabel('Percentate of Variance Explained')
    plt.xlabel('Principal Component')
    plt.title ('PCA Scree Plot')
    plt.show()

makePCA(tissueMatrix_norm)







#validator.validate(phiList,prediction,logScale="yes",xLabel="Phi",yLabel="Prediction")




