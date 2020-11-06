#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:10:39 2020

@author: lu
"""

#Napping genes with same gene name

fn1="/home/lu/AcrossTissue/csvs/fiveStageMeanEmpiricalExp.csv"
fn2="/home/lu/AcrossTissue/csvs/wholeGenomePhi.csv"

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:50:09 2019

@author: lu
"""



def getExpDictFromFile(fname): 
    f1=open(fname,"r+")
    lines=f1.readlines()
    geneExpressionMap=dict()
    #col 0 is the gene id, then the corespoding expression value at col1, then col 1 is gene id, and expression value at col2 and so on...
    for line in lines:
        splitList=line.rstrip().split(",")
        try:
            geneId=splitList[0].rstrip().lower()
            expression=splitList[1]
    #        print "geneId: %s"%geneId
    #        print "Expression: %s"%expression
        except:
            print "index 1 no expression"
            continue
        if "null" in geneId:
            continue
        if expression=='null' or expression=="NA": # here, I treat the gene with null as expression measurement to have a expression of 0.01
    #                    expression=2
            continue #skip those with null as expression
        #skip the first one
        if geneId in geneExpressionMap:
            continue
        else:
            geneExpressionMap[geneId]=(expression)
                
    return geneExpressionMap

geneExpDict=getExpDictFromFile(fn1)
phiDict=getExpDictFromFile(fn2)


geneNameList=[]
expList=[]
phiList=[]


for geneName in geneExpDict:
    if geneName in phiDict:
        geneNameList.append(geneName)
        expList.append(float(geneExpDict[geneName]))
        phiList.append(float(phiDict[geneName]))
      
import validator

validator.validate(phiList,expList,logScale="no",xLabel="phi",yLabel="EmpiricalExpression")   
cnt=0
for item in expList:
    if item==0:
        cnt+=1
print cnt
        
