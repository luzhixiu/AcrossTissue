#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 02:09:08 2020

@author: lu
"""
import statistics as ss
import readCSVfiles as RCF
import numpy as np
import matplotlib.pyplot as plt
from math import log
import findSequenceById

topCutOff=3000
fname="../csvs/LifeStageGenes_collapse.csv"


raw_matrix,rowHeaderList,columnHeaderList= RCF.readCSV(fname,readAsString=True)


raw_matrix=np.concatenate(([columnHeaderList],raw_matrix),axis=0)
raw_matrix_transposed=np.transpose(raw_matrix)


sortList=raw_matrix_transposed[:,5]

sortList=[float(x) for x in sortList]

sortIndex=np.asarray(sortList).argsort()
sortIndex=np.flip(sortIndex)

matrix=raw_matrix_transposed[sortIndex]

matrix=(matrix[:topCutOff])


sequenceDict=findSequenceById.findSequenceByID("../Fastas/c_elegan.fasta",idType="Gn")



LS_WBID_Dict=dict()

for lst in matrix:
    WBID=lst[0]
    LS=lst[1]
    if LS in LS_WBID_Dict:
        LS_WBID_Dict[LS].append(WBID)
    else:
        LS_WBID_Dict[LS]=[WBID]
print(LS_WBID_Dict.keys())


for LS in LS_WBID_Dict:
    prefix="../Fastas/"
    fileName=LS.replace(" ","_")
    f=open(prefix+fileName+".fasta","w")
    WBIDList=LS_WBID_Dict[LS]
    print(len(WBIDList))

    for WBID in WBIDList:
        if WBID in sequenceDict:
            sequence=sequenceDict[WBID]
            f.write(">"+WBID+"\n")
            f.write(sequence+"\n")
    f.close()
            
    
        


