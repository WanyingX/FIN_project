#!/bin/python

# This script is used to integrate huge loop data from two different conditions.

import os
import sys
import pandas as pd
import numpy as np

def readfile(filename):
	df=pd.read_csv(filename,sep='\t',names=['a1','a2','ratio'],header=None,index_col=False)
#	df['a1'] = df['a1'].map(lambda x: x.replace('A_', ''))
#	df['a2'] = df['a2'].map(lambda x: x.replace('A_', ''))
	df.index=df['a1']+':'+df['a2']	
	return df


name1=sys.argv[1]
name2=sys.argv[2]
name3=sys.argv[3]

df1=readfile(name1+".merged")
df2=readfile(name2+".merged")
#gt1=readfile(name3+".lifted")
gt1=pd.read_csv(name3,sep='\t',names=['a1','a2'],header=None,index_col=False)
#gt2=readfile("lostloops.lifted")
#gt3=readfile('commonloops.lifted')
index1=gt1['a1']+':'+gt1['a2']
#index1.extend(gt2.index.to_list())
#index1=set(index1)

f1=df1[:][df1.index.isin(index1)]
f2=df2[:][df2.index.isin(index1)]

#f1.to_csv("gainloops.lifted.wt",sep='\t',index=False)

f1['KO']=f2['ratio']
f1.to_csv(name3+".lifted.ratio",sep='\t',index=False)

