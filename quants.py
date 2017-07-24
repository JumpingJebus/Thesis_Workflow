#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 20:00:02 2017

@author: stephen
"""

import pandas as pd

pd.options.mode.chained_assignment = None

for t in range(51,67):
    quants= pd.read_table("/home/stephen/Downloads/quants/salmon/SRR33740" + str(t) + "/quant.sf")
    for i in range(0,(len(quants["Name"]))):
        quants["Name"][i] = quants["Name"][i][3:20]
    quants.to_csv("/home/stephen/Downloads/quants/salmon/SRR33740" + str(t) + "/quants1.sf", sep="\t", encoding = "utf-8",index= False)
