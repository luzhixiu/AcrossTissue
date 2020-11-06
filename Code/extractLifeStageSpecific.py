#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:27:27 2020

@author: lu
"""

f=open("/home/lu/Desktop/development_ontology.WS274.obo","r+")
lines=f.readlines()
f.close()
f=open("wormbase_development_extract.txt","w+")
for  i in range(len(lines)):
    line=lines[i]
    if "[Term]" in line:
        wb_id= lines[i+1]
        lifestage= lines[i+2]
        print wb_id
        print lifestage
        f.write (wb_id)
        f.write (lifestage+"\n")

