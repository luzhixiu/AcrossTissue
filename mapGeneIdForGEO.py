#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:50:09 2019

@author: lu
"""

f=open("caiAndCAI.csv","r+")
lines=f.readlines()

header=lines[0]
print header
lines=lines[1:]


geneExpressionMap=dict()
cnt=0
#col 0 is the gene id, then the corespoding expression value at col1, then col 1 is gene id, and expression value at col2 and so on...
for line in lines:
    splitList=line.rstrip().split(",")
    print splitList
    for i in range(len(splitList)):
            if i%2==0:
                try:
                    geneId=splitList[i].rstrip()
                    expression=splitList[i+1]
                except:
                    print "index %d no expression"%i
                if "null" in geneId:
                    continue
                if expression=='null' or expression=="NA": # here, I treat the gene with null as expression measurement to have a expression of 0.01
#                    expression=2
                    continue #skip those with null as expression
                if geneId in geneExpressionMap:
                    geneExpressionMap[geneId].append(expression)
                else:
                    geneExpressionMap[geneId]=[]
                    geneExpressionMap[geneId].append(expression)
                    
                 

        
f=open("cai_comparison.csv","w+")            
for key in geneExpressionMap:
    if not len(geneExpressionMap[key])==2:
        continue
    f.write(key)
    f.write(",")    
    for expression in geneExpressionMap[key]:
        f.write(str(expression))
        f.write(",")
    f.write("\n")
f.close()