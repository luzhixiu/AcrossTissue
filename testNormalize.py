#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:15:10 2020

@author: lu
"""
from scipy.stats.mstats import gmean 
from sklearn.preprocessing import RobustScaler

testList=[1,1,1,9,9,9]




def scaleToOne(lst):
    avg=gmean(lst)
    result=[x/avg for x in lst]
    print gmean(result)
