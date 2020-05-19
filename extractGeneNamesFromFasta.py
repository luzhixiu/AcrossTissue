#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 13:01:47 2020

@author: lu
"""

def findSequenceByID(inputFile,idType="locus_tag"):
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

geneNameList=findSequenceByID("/home/lu/Desktop/sequences/c_elegan.fasta",idType="gene")
print len(geneNameList)



import csv

with open('geneNames.csv','wb') as result_file:
    wr = csv.writer(result_file, dialect='excel')
    for item in geneNameList:
        wr.writerow([item,])