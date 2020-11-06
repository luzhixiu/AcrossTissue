#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 12:09:29 2020

@author: lu
"""
import pandas as pd
import seaborn as sns;

df = pd.read_csv("/home/lu/AcrossTissue/csvs/cEl.csv")


sns.lmplot(x='Attack', y='Defense', data=df)