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

rd.seed(0)
tissueSamples=25


sequenceDict=findSequenceById.findSequenceByID("/home/lu/Desktop/sequences/c_elegan.fasta",idType="gene")
geneNameList=[]

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

f=open("cEl.csv","r+")
lines=f.readlines()

#21 tissues or lifestages

tissueMatrix=[]
phiList=[]

for i in range(tissueSamples):
    tissueMatrix.append([])

for line in lines:
    
#..............This line is for fasting testings, uncomment when running    
#    randomN=rd.randrange(1,101)
#    if randomN>=5:
#        continue
    line=line.rstrip()
    splitList=line.split(",")
    geneName=splitList[0]
    if geneName in sequenceDict:
        geneNameList.append(geneName)
        sequence=sequenceDict[geneName]
        codonList=loadSequence(sequence)
    else:
        continue
    for i in range(1,len(splitList)):
        if len(splitList[i])==0:
            tissueMatrix[i-1].append(0)
        else:
            tissueMatrix[i-1].append(float(splitList[i]))

tissueMatrix=np.asarray(tissueMatrix)
indexList=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

tissueMatrix=tissueMatrix[indexList]
print tissueMatrix.shape

print len(tissueMatrix)



codontable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

codonList=codontable.keys()








def corelate(list1,list2):
    return ss.pearsonr(list1,list2)[0]


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
                  

def average(lst):
    return float(sum(lst))/len(lst)


def scaleToOne(lst):
    avg=average(lst)
    result=[x/avg for x in lst]
    return result

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

sampleLables=["hermaphrodite, NSM, L1 larva Ce","hermaphrodite, gonad, adult Ce", "hermaphrodite, motor neuron, L2 larva Ce", "hermaphrodite, neuron, L1 larva Ce" ,"hermaphrodite, organism, 3-fold embryo Ce","hermaphrodite, organism, 4-cell embryo Ce","hermaphrodite, organism, L1 larva Ce","hermaphrodite, organism, L2 larva Ce","hermaphrodite, organism, L2d-dauer molt Ce","hermaphrodite, organism, L3 larva Ce","hermaphrodite, organism, L4 larva Ce","hermaphrodite, organism, adult Ce","hermaphrodite, organism, dauer larva Ce","hermaphrodite, organism, elongating embryo Ce","hermaphrodite, organism, enclosing embryo Ce","hermaphrodite, organism, fully-elongated embryo Ce","hermaphrodite, organism, gastrulating embryo Ce","hermaphrodite, organism, late cleavage stage embryo Ce","hermaphrodite, organism, newly molted young adult hermaphrodite Ce","hermaphrodite, organism, post dauer stage Ce","hermaphrodite, organism, proliferating embryo Ce","hermaphrodite, pharyngeal muscle cell, fully-elongated embryo Ce","hermaphrodite, somatic cell, embryo Ce","male, organism, L4 larva Ce","male, organism, embryo Ce"
]
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



#find the percent of a cutoff in a list, returns the value at that cutoff,
#f should be something like 0.25, which is the cut off for bot 25 percent or top 75
def findPercent(lst,f):
    lst=sorted(lst)
    return lst[int(len(lst)*f)]

def logify(lst):
    logList=[]
    for i in lst:
        if i==0:
            i=0.00001
        logList.append(np.log(i))
    return logList
    

#print phiList

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
    

tissueMatrix_norm=[]
for tissueList in tissueMatrix:
    tissueMatrix_norm.append(scaleToOne(tissueList))
    
tissueMatrix_log_norm=[]
for tissueList in tissueMatrix_norm:
    tissueMatrix_log_norm.append(logify(tissueList))
    

tissueMatrix_norm_trans=np.transpose(tissueMatrix_norm)
chosenLables=np.asarray(sampleLables)[indexList]


genomeList=[]




for k in range(len(tissueMatrix_log_norm)):
    wholeGenomeSeq=""
    tissueList=tissueMatrix_log_norm[k]
#for tissueList in tissueMatrix_norm:
    geneNameSelect=[]
    cutOff=findPercent(tissueList,0.9)
    for i in range(len(tissueList)):
        if tissueList[i]>=cutOff:
            wholeGenomeSeq+= sequenceDict[geneNameList[i]]
    genomeList.append(wholeGenomeSeq)
print len(genomeList)
    



import calculateRSCU as RSCU

SynonymousCodons ={
'CYS': ['TGT', 'TGC'],
'ASP': ['GAT', 'GAC'],
'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
'GLN': ['CAA', 'CAG'],
'MET': ['ATG'],
'ASN': ['AAC', 'AAT'],
'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
'LYS': ['AAG', 'AAA'],
'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
'PHE': ['TTT', 'TTC'],
'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
'ILE': ['ATC', 'ATA', 'ATT'],
'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 'HIS': ['CAT', 'CAC'],
'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
'TRP': ['TGG'],
'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
'GLU': ['GAG', 'GAA'],
'TYR': ['TAT', 'TAC']}

#f=open("CUB_ByLifeStage.txt","w")
#for i in range(len(genomeList)):
#    genomeString=genomeList[i]
#    result=RSCU.getRSCU(genomeString)
#    f.write(chosenLables[i])
#    f.write("\n")
#    for j in result:
#        for t in j:
#            print str(t)
#            f.write(str(t))
#            f.write("\n")
#f.close()
##    for key in SynonymousCodons:
##        print key
#    


emb=tissueMatrix[4]
larvae_L1=tissueMatrix[6]
L3_L4=tissueMatrix[9]
dauer=tissueMatrix[12]
adult=tissueMatrix[11]


meanExpressionList=[]


embList=[]
larvaList=[]
L3_L4List=[]
dauerList=[]
adultList=[]
allList=[]

embExp=[]
larvaeExp=[]
L3L4Exp=[]
dauerExp=[]
adultExp=[]


print len(emb)
print len(geneNameList)
for i in range(len(emb)):
    allList.append(geneNameList[i])
    mean=average([emb[i],larvae_L1[i],L3_L4[i],dauer[i],adult[i]])
    threshold=2*mean
    if emb[i]>threshold:
        embList.append(geneNameList[i])
        embExp.append(emb[i])
    if larvae_L1[i]>threshold:
        larvaList.append(geneNameList[i])    
        larvaeExp.append(larva_L1[i])
    if L3_L4[i]>threshold:
        L3_L4List.append(geneNameList[i])  
        L3L4Exp.append(L3_L4[i])
    if dauer[i]>threshold:
        dauerList.append(geneNameList[i])
        dauerExp.append(dauer[i])                
    if adult[i]>threshold:
        adultList.append(geneNameList[i])            
        adultExp.append(adult[i])
#sequenceDict=findSequenceById.findSequenceByID("/home/lu/Desktop/sequences/c_elegan.fasta",idType="gene")

def writeToFasta(geneTagList,fName):
    global sequenceDict
    f=open(fName,"w+")
    for i in range(len(geneTagList)):
        f.write(">%d %s\n"%(i,geneTagList[i]))
#        sequence=sequenceDict[geneTagList[i]]
#        f.write(sequence)
#        f.write("\n")
        
#writeToFasta(embList,"emb.fasta")
#writeToFasta(larvaList,"larvae.fasta")
#writeToFasta(L3_L4List,"L3_L4.fasta")
#writeToFasta(dauerList,"dauer.fasta")
#writeToFasta(adultList,"adult.fasta")
writeToFasta(allList,"wholeGenome.seqName")





