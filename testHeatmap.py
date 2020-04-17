#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 10:01:18 2020

@author: lu
"""

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
uniform_data = np.random.rand(10000, 20)
print uniform_data.shape
ax = sns.heatmap(uniform_data)
