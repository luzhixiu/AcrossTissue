#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 11:53:50 2020

@author: lu
"""
import numpy as np
from sklearn.decomposition import PCA


def makePCA(X,n=3):
    

    pca = PCA(n_components=2)
    pca.fit(X)
    PCA(n_components=2)
    print pca.explained_variance_ratio_