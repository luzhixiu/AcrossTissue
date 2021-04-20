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
        if started:
            actualCodonList.append(codon)
        if codon==startCodon:
            started=True
    codonList=actualCodonList
   # print "codon readed successful, the number of codon in this sequence is %d"%(len(codonList))
    return codonList


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
    
# key is the aa sequence and value is a list of codon sequences that produce same aa
sameAAdiffCodon=dict()
for geneName in sequenceDict:
    seq=sequenceDict[geneName]
    protSeq=geneSequenceToProteinSequence(seq)
    if protSeq in sameAAdiffCodon:
        if seq not in sameAAdiffCodon[protSeq]:
            sameAAdiffCodon[protSeq].append(seq)
    else:
        sameAAdiffCodon[protSeq]=[seq]

for protSeq in sameAAdiffCodon:
    codonCandidateList=sameAAdiffCodon[protSeq]
    if len(codonCandidateList)>=2:
        print ("AA Sequence: ")
        print (protSeq)
        print ("Codon Sequence: ")
        print (codonCandidateList)
        
    

    














