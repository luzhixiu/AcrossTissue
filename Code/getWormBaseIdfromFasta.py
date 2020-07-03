#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 23:24:07 2020

@author: lu
"""

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


wormbaseNames=findSequenceByID("../Fastas/c_elegan.fasta",idType="Gn")
import csv

with open('wormBaseId.csv','wb') as result_file:
    wr = csv.writer(result_file, dialect='excel')
    for item in wormbaseNames:
        wr.writerow([item,])

#print wormbaseNames
