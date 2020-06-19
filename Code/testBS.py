#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:48:15 2020

@author: lu
"""

lst=[3,7,13,67,322]


def testBS(lst,low,high,target):
    print low
    if low==high:
        return lst[low]
    if (low+high)%2==0:
        mid=(low+high)/2
    else:
        mid=(low+high-1)/2
    if lst[mid]>=target:
        testBS(lst,low,mid,target)
    else:
        testBS(lst,mid,high,target)
    
        
    
print testBS(lst,0,len(lst)-1,7)