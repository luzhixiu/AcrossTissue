#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:10:24 2020

@author: lu
"""

import matplotlib.pyplot as plt
import Bio.SeqUtils.CodonUsage  as BIOCUB
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
from sklearn.preprocessing import minmax_scale
import findSequenceById
import numpy as np
import scipy.stats as ss
import math
import random as rd
import statsmodels.api as sm
from Bio import SeqIO

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
        if codon==startCodon:
            if started:
                break
            else:
                started=True        
        if started:
            actualCodonList.append(codon)

   # print "codon readed successful, the number of codon in this sequence is %d"%(len(codonList))
    return actualCodonList


sequenceDict=findSequenceById.findSequenceByID("/home/lu/AcrossTissue/Fastas/c_elegan.fasta",idType="gene")



CodonsDict = {
    "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
    "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
    "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
    "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
    "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
    "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
    "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
    "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
    "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
    "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
    "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
    "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
    "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
    "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
    "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
    "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}

    
inverseTable=  {
    "CYS": ["TGT", "TGC"],
    "ASP": ["GAT", "GAC"],
    "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "GLN": ["CAA", "CAG"],
    "MET": ["ATG"],
    "ASN": ["AAC", "AAT"],
    "PRO": ["CCT", "CCG", "CCA", "CCC"],
    "LYS": ["AAG", "AAA"],
    "STOP": ["TAG", "TGA", "TAA"],
    "THR": ["ACC", "ACA", "ACG", "ACT"],
    "PHE": ["TTT", "TTC"],
    "ALA": ["GCA", "GCC", "GCG", "GCT"],
    "GLY": ["GGT", "GGG", "GGA", "GGC"],
    "ILE": ["ATC", "ATA", "ATT"],
    "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "HIS": ["CAT", "CAC"],
    "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "TRP": ["TGG"],
    "VAL": ["GTA", "GTC", "GTG", "GTT"],
    "GLU": ["GAG", "GAA"],
    "TYR": ["TAT", "TAC"],
}

codonTable=dict()
for aa in inverseTable:
    codonList=inverseTable[aa]
    for codon in codonList:
        codonTable[codon]=aa
    
    


def geneSequenceToProteinSequence(seq):
    codonList=loadSequence(seq)
    aaString=""
    for codon in codonList:
        aa=codonTable[codon]
        aaString+=aa
    return aaString
    
# # key is the aa sequence and value is a list of codon sequences that produce same aa
# sameAAdiffCodon=dict()
# for geneName in sequenceDict:
#     seq=sequenceDict[geneName]
#     protSeq=geneSequenceToProteinSequence(seq)
#     if protSeq in sameAAdiffCodon:
#         if seq not in sameAAdiffCodon[protSeq]:
#             sameAAdiffCodon[protSeq].append(seq)
#     else:
#         sameAAdiffCodon[protSeq]=[seq]

# for protSeq in sameAAdiffCodon:
    
#     codonCandidateSeqList=sameAAdiffCodon[protSeq]
#     if len(codonCandidateSeqList)>=2:
#         indexDict=dict()
#         codonCandidateList=[loadSequence(seq) for seq in codonCandidateSeqList] 
#         for i in range(len(codonCandidateList[0])):
#             codonOptions=[]
#             for codonList in codonCandidateList:
#                 codon=codonList[i]
#                 codonOptions.append(codon)
#             if len(set(codonOptions))>1:
#                 indexDict[i]=codonOptions
#         print(indexDict)
                
#n being the number of codons, 2 for codon pair, 3 for codon triplets and so on
def codonListToPairList(codonList,n):
    kCodonList=[]
    #this part can be optimized if it gets slow* 
    for i in range(0,len(codonList)-n+1):
        kCodon=""
        for j in range(i,i+n):
            kCodon+=codonList[j]
        if len(kCodon)==n*3:
            kCodonList.append(kCodon)
    return kCodonList
    

codonPairList_All=[]
for geneName in sequenceDict:
    seq=sequenceDict[geneName]
    codonList=loadSequence(seq)
    codonPairList=codonListToPairList(codonList,2)
    codonPairList_All.append(codonPairList)
    
codonPairList_Flat= [item for sublist in codonPairList_All for item in sublist]
            
from collections import Counter
import operator
codonPairCntDict_celegan=Counter(codonPairList_Flat)
codonPairCntDict_sorted = dict( sorted(codonPairCntDict_celegan.items(), key=operator.itemgetter(1),reverse=True))

codonKeys=codonPairCntDict_sorted.keys()

codonPairCntDict_sorted_celegan=codonPairCntDict_sorted 

        
    

    

sequenceDict=findSequenceById.findSequenceByID("/home/lu/AcrossTissue/Fastas/s288c.fasta")



codonPairList_All=[]
for geneName in sequenceDict:
    seq=sequenceDict[geneName]
    codonList=loadSequence(seq)
    codonPairList=codonListToPairList(codonList,2)
    codonPairList_All.append(codonPairList)
    
codonPairList_Flat= [item for sublist in codonPairList_All for item in sublist]
            
from collections import Counter
import operator
codonPairCntDict_yeast=Counter(codonPairList_Flat)

codonPairCntDict_sorted_yeast = dict( sorted(codonPairCntDict_yeast.items(), key=lambda x: x[0].lower()) )

codonPairCntDict_sorted_cel = dict( sorted(codonPairCntDict_celegan.items(), key=lambda x: x[0].lower()) )

celList=[]
yeastList=[]



for key in codonPairCntDict_yeast:
    yeastList.append(codonPairCntDict_yeast[key])
    celList.append(codonPairCntDict_celegan[key])

import scipy.stats as ss



def testCorelation(x,y,corelationFunction):

    if "pearson" in corelationFunction:
#        print ("p value %f"%stats.pearsonr(x, y)[1])
        return  ss.pearsonr(x, y)[0]
    elif "spearman" in corelationFunction:
        return ss.spearmanr(x,y,nan_policy="omit")[0]
    elif "kendall" in corelationFunction:
        return ss.kendalltau(x,y,nan_policy="omit")[0]
    
spearmanr=testCorelation(yeastList,celList,"spearman")
    
import validator

validator.plotCorelation(yeastList, celList,xLabel="Yeast",yLabel="C. Elegan")

yeastKeyList=[]
celKeyList=[]
for key in codonPairCntDict_yeast:
    yeastList.append(codonPairCntDict_yeast[key])
    celList.append(codonPairCntDict_celegan[key])


codonPairCntDict_celegan_sorted_value=dict(sorted(codonPairCntDict_celegan.items(), key=lambda item: item[1],reverse=True))
codonPairCntDict_yeast_sorted_value=dict(sorted(codonPairCntDict_yeast.items(), key=lambda item: item[1],reverse=True))

celegan_ranked_codons=list(codonPairCntDict_celegan_sorted_value.keys())
yeast_ranked_codons=list(codonPairCntDict_yeast_sorted_value.keys())

# print(testCorelation(celegan_ranked_codons,yeast_ranked_codons,"spearmanr"))

celegan_codon_pairt_rank_dict={k: v for v, k in enumerate(celegan_ranked_codons)}
yeast_codon_pairt_rank_dict={k: v for v, k in enumerate(yeast_ranked_codons)}

rankDiffDict=dict()
for key in celegan_codon_pairt_rank_dict:
    rankDiff=celegan_codon_pairt_rank_dict[key]-yeast_codon_pairt_rank_dict[key]
    rankDiffDict[key]=rankDiff

diffList=list(rankDiffDict.values())
plt.hist(diffList)
    