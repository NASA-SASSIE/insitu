#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 19:14:06 2022

@author: suzanne
"""

import pandas as pd
import numpy as np

outHeaders = {'Date':'Date',
           'DateTimeStr':'Date',
           'CTD-S2':'Salinity',
           'S1':'Salinity',
           'Salinity-0':'Salinity',
           'Ts':'Temperature',
           'T1':'Temperature',
           'WaterTemp-0':'Temperature'}
print(outHeaders['S1'])
exit(-1)

num=641
print(bin(num)[2:].zfill(16)[::-1])
exit(-1)
binnum = bin(num)[2:].zfill(16)[::-1]
print(binnum)
print()
print(binnum[::-1])
exit(-1)
if binnum[5]:
    print('roughness correction')
if binnum[7]:
    print('land')
if binnum[8]:
    print('ice')
exit(-1)

# df1 = pd.DataFrame({'col1':[11,12], 'col2':[23,24]})
# print(df1.head())
      
# df2 = pd.DataFrame({'colA':['ant','aardvark'], 
#                     'colB':['bat','badger'],
#                     'colC':['cat','cougar']
#                     })
# print(df2.head())
      
# filecsv = 'testfile.xls'

# with pd.ExcelWriter(filecsv) as writer:
#     df1.to_excel(writer, sheet_name='UTO 1',index=False)
#     df2.to_excel(writer, sheet_name='Swift 12',index=False)
    
      
      
      