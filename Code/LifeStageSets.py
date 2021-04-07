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


# only one of these preset options should be specified, the "foldDiffCutOff"
# simply changes the cutoff and filter out the genes below this number, 
# "topFromEachGroup" select the same number of genes from each group (LS).
foldDiffCutOff=0
topFromEachGroup=300
####==========================================
inputFileName="../csvs/LSB.csv"




def main():
    global foldDiffCutOff,topFromEachGroup
    global inputFileName
    prefix="../Fastas/"
    
    
    raw_matrix,rowHeaderList,columnHeaderList= RCF.readCSV(inputFileName,readAsString=True)
    raw_matrix=np.concatenate(([columnHeaderList],raw_matrix),axis=0)
    raw_matrix_transposed=np.transpose(raw_matrix)


    if foldDiffCutOff>0:
        postfix_name="_foldDiff_"+str(foldDiffCutOff)
        print("Using foldDiff Cut Off of ",foldDiffCutOff)
        selectedGeneNumbers=0

        foldDiffList=raw_matrix_transposed[:,5]
        foldDiffList=[float(x) for x in foldDiffList]
        sortIndex=np.asarray(foldDiffList).argsort()
        sortIndex=np.flip(sortIndex)
        print(sortIndex)
        selectedGeneNumbers=sorted(foldDiffList,reverse=True).index(foldDiffCutOff)
        print("Number of Genes Selected: ",selectedGeneNumbers)
        matrix=raw_matrix_transposed[sortIndex]
        matrix=(matrix[:selectedGeneNumbers])
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
            fileName=LS.replace(" ","_")
            f=open(prefix+fileName+postfix_name+".fasta","w")
            WBIDList=LS_WBID_Dict[LS]
            print(len(WBIDList))
            cnt=0
            for WBID in WBIDList:
                if WBID in sequenceDict:
                    sequence=sequenceDict[WBID]
                    f.write(">"+WBID+"\n")
                    f.write(sequence+"\n")
                else:
                    cnt+=1
            print("Numbers of genes that are not found in genome: ", cnt)
            f.close()
    else:
        print("Using The Top",topFromEachGroup,"Genes From Each Group")
        postfix_name="_TopFromEachGroup_"+str(topFromEachGroup)
        foldDiffList=raw_matrix_transposed[:,5]
        foldDiffList=[float(x) for x in foldDiffList]
        sortIndex=np.asarray(foldDiffList).argsort()
        sortIndex=np.flip(sortIndex)
        matrix=raw_matrix_transposed[sortIndex]
        
        sequenceDict=findSequenceById.findSequenceByID("../Fastas/c_elegan.fasta",idType="Gn")
        LS_WBID_Dict=dict()
        LS_Set=set()
        for lst in matrix:
            LS=lst[1]
            LS_Set.add(LS)
        for ls in LS_Set:
            print(ls)
            fileName=ls.replace(" ","_")
            f=open(prefix+fileName+postfix_name+".fasta","w")
            cnt=0
            for lst in matrix:
                WBID=lst[0]
                LS=lst[1]
                if LS==ls:
                    if cnt<topFromEachGroup:
                        if WBID in sequenceDict:
                            cnt+=1
                            sequence=sequenceDict[WBID]
                            f.write(">"+WBID+"\n")
                            f.write(sequence+"\n")

                    else:
                        break
            f.close()                           
                    
                
            
    
        
main()

