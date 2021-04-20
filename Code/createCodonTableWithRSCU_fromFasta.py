#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 23:02:08 2021

@author: lu
"""


import matplotlib.pyplot as plt
import Bio.SeqUtils.CodonUsage  as BIOCUB
import numpy as np
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
import findSequenceById
from Bio import SeqIO


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



class customizedCUB(CodonAdaptationIndex):
   
    def _count_codons(self, fasta_file):
        cnt=0
        with open(fasta_file) as handle:

            # make the codon dictionary local
            self.codon_count = CodonsDict.copy()

            # iterate over sequence and count all the codons in the FastaFile.
            for cur_record in SeqIO.parse(handle, "fasta"):
                goodSequence=True
                # make sure the sequence is lower case
                if str(cur_record.seq).islower():
                    dna_sequence = str(cur_record.seq).upper()
                else:
                    dna_sequence = str(cur_record.seq)

                for i in range(0, len(dna_sequence), 3):
                    codon = dna_sequence[i : i + 3]
                    if codon not in self.codon_count:
                        goodSequence=False
                
                if goodSequence:
                    cnt+=1
                    for i in range(0, len(dna_sequence), 3):
                        codon = dna_sequence[i : i + 3]
                        if codon in self.codon_count:
                            self.codon_count[codon] += 1
        print("Number of Sequences used to create RSCU: ",cnt)
        
                        
CUB=customizedCUB()   
        

CUB.generate_index("/home/lu/AcrossTissue/Fastas/random_c_elegan_300.fasta")

RSCU_Dict=(CUB.index)


figDim=5
fig, axs = plt.subplots(figDim, figDim,figsize=(15,15))
        
c=0
r=0        


for aa in inverseTable:
    colorList=list('rgbkymc')
    print(c,r)
    codonList=inverseTable[aa]
    
    RSCU_List=[RSCU_Dict[x] for x in codonList]
    sort_index = np.flip(np.argsort(RSCU_List))
    codonList=[codonList[i] for i in sort_index]        
    RSCU_List=[RSCU_List[i] for i in sort_index]
    colorList=[colorList[i] for i in sort_index]
    axs[r,c].bar(codonList,RSCU_List,color=colorList)
    axs[r,c].set_title(aa)


    c+=1
    if c>=figDim:
        r+=1
        c=0
    
fig.tight_layout()      
    
    
    
