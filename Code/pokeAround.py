#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:10:24 2020

@author: lu
"""

import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale
import findSequenceById
import numpy as np
import scipy.stats as ss
import math
import random as rd
import statsmodels.api as sm

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


sequenceDict=findSequenceById.findSequenceByID("/home/lu/Desktop/sequences/c_elegan.fasta",idType="gene")


    
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

inverseTable=  {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'], 
                'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
              'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'], 
              'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'], 
              'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'], 
              'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'], 
              'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}    
    


def geneSequenceToProteinSequence(seq):
    codonList=loadSequence(seq)
    aaString=""
    for codon in codonList:
        aa=codontable[codon]
        aaString+=aa
    return aaString
    
#key is the aa sequence and value is a list of codon sequences that produce same aa
sameAAdiffCodon=dict()
for geneName in sequenceDict:
    seq=sequenceDict[geneName]
    protSeq=geneSequenceToProteinSequence(seq)
    if protSeq in sameAAdiffCodon:
        sameAAdiffCodon[protSeq].append(seq)
    else:
        sameAAdiffCodon[protSeq]=[seq]

for protSeq in sameAAdiffCodon:
    codonCandidateList=sameAAdiffCodon[protSeq]
    if len(codonCandidateList)>10:
        print "AA Sequence: "
        print protSeq
        print "Codon Sequence: "
#        print len(codonCandidateList)
        for codonSeq in codonCandidateList:
            print codonSeq


    


    
    

sequence="ATGAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTTAGTAATGA"
#print geneSequenceToProteinSequence(sequence)