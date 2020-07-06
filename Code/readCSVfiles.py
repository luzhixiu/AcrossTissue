#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 13:39:22 2020

@author: lu
"""


import numpy as np
#read in a csv files into matrix, default assumes there is a rowHeader and return it, 
#if hasrowHeader is passed in as false, return empty rowHeaderlist
#input: csv file path
#output: matrix of the csx, each column is a list relative to coresponding rowHeader location
#for example, matrix[3] is a list from responding column, and rowHeader[3] describes this column 
def readCSV(fname, hasrowHeader=True,hasColHeader=True,partialIndexList=[]):
    rowHeader=[]
    columnHeaderList=[]
    matrix=[]
    f=open(fname,"r+")
    lines=f.readlines()
    rowHeader=lines[0]
    lines=lines[1:]
    rowHeaderList=rowHeader.strip().split(",")
    nSamples=len(rowHeaderList)
    if hasColHeader:
        colSize=nSamples-1
    else:
        colSize=nSamples
    for i in range(colSize):
        matrix.append([])
    
    for line in lines:
        line =line.strip()
        if len(line)==0:
            continue
    #This line is for faster testings, comment when running    
    #    randomN=rd.randrange(1,101)
    #    if randomN>=5:
    #        continue
        splitList=line.split(",")
        if len(splitList)==0:
            continue
        
        columnHeaderList.append(splitList[0])
        for i in range(1,len(splitList)):
            if len(splitList[i])==0:
                matrix[i-1].append(0)
            else:
                matrix[i-1].append(float(splitList[i]))
    #if the partialIndexList is given, only return partial matrix
    if not len(partialIndexList) == 0:        
        matrix=[matrix[index] for index in partialIndexList]
        rowHeaderList=[rowHeaderList[index] for index in partialIndexList]
        columnHeaderList=[columnHeaderList[index] for index in partialIndexList]
    return matrix,rowHeaderList,columnHeaderList
        
        
    
#result=readCSV("/home/lu/AcrossTissue/csvs/cEl.csv")
#print result