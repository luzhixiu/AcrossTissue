#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 13:39:22 2020

@author: lu
"""



#read in a csv files into matrix, default assumes there is a header and return it, 
#if hasHeader is passed in as false, return empty headerlist
#input: csv file path
#output: matrix of the csx, each column is a list relative to coresponding header location
#for example, matrix[3] is a list from responding column, and header[3] describes this column 
def readCSV(fname, hasHeader=True):
    header=[]
    matrix=[]
    
    f=open(fname,"r+")
    lines=f.readlines()
    header=lines[0]
    headerList=header.rstrip.split(",")
    nSamples=len(headerList)
    
    for i in range(nSamples):
        matrix.append([])
    
    for line in lines:
    #This line is for faster testings, comment when running    
    #    randomN=rd.randrange(1,101)
    #    if randomN>=5:
    #        continue
        splitList=line.rstrip().split(",")
        for i in range(1,len(splitList)):
            if len(splitList[i])==0:
                matrix[i-1].append(0)
            else:
                matrix[i-1].append(float(splitList[i]))
    return matrix,header
        
        
    
result=readCSV("/home/lu/AcrossTissue/csvs/cEl.csv")
print result