#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 23:24:07 2020

@author: lu
"""
import csv
import string
def findSequenceByID(inputFile,idType="locus_tag"):
    print ("Selected id Type: %s"%(idType))
    wormbaseIdList=[]
    from Bio import SeqIO
    records=SeqIO.parse(inputFile, "fasta")
    
    cnt=0
    mySum=0
    for record in records:
        mySum+=1
        header=str(record.description)
        startTargetIndex=header.find(str(idType))
        
        if startTargetIndex<0:
#            print "couldn't find the target idType"
            cnt+=1
            wormbaseIdList.append("NA")
            continue
        startIndex=startTargetIndex+len(idType)+1
        idName=""
        charIndex=startIndex
        while not (header[charIndex]=="]" or header[charIndex]==","):
            idName+=header[charIndex]
            charIndex+=1
        wormbaseIdList.append(idName)

    return wormbaseIdList

def findSequenceByID2(inputFile,idType="locus_tag"):
    geneNameList=[]
    geneDict=dict()
    from Bio import SeqIO
    records=SeqIO.parse(inputFile, "fasta")
    cnt=0
    mySum=0
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
#        if idName not in geneDict:
#            geneDict[idName]=str(record.seq)
        geneNameList.append(idName)
        
    return geneNameList

wormbaseNames=findSequenceByID("../Fastas/c_elegan.fasta",idType="Gn")
geneNameList=findSequenceByID2("../Fastas/c_elegan.fasta",idType="gene")





#
with open('wormBaseId.csv','w') as result_file:
    wr = csv.writer(result_file, dialect='excel')
    for i in range(len(wormbaseNames)):
        wr.writerow([wormbaseNames[i],geneNameList[i]])
        
idMap=dict()
for i  in range(len(wormbaseNames)):
    idMap[wormbaseNames[i]]=(geneNameList[i])
print (idMap)

import readCSVfiles as RCF

matrix,header1,header2=RCF.readCSV("../csvs/6LS_collapse.csv")


geneNameList=[]
for WBID in header2:
    
    if WBID in idMap:
        geneName=idMap[WBID]
        geneNameList.append(geneName)
    else: 
        print("No found: ",geneName)
        geneNameList.append("NA")

print(geneNameList)
    
import csv

with open('temp.csv','w') as result_file:
    wr = csv.writer(result_file, dialect='excel')
    for item in geneNameList:
        wr.writerow([item,])    
    
    
    
    
    
    
    
    