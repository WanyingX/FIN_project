#!/bin/python
# This script is used to lift over chromatin loop from different Hi-C experiments.

import pandas as pd
import sys
import os


bed={}

with open(sys.argv[1],'r') as f:

	for line in f:
		f1, a1 = line.rstrip().split('\t')
		if f1 not in bed:
			bed[f1]=a1
f.close()

dic={}

for line in sys.stdin:
	f1, f2 = line.rstrip().split('\t')[0:2]
	if f1 in bed and f2 in bed:
		a1, a2 = bed[f1],bed[f2]
		if a1 not in dic:
			dic[a1]={}

		if a2 not in dic[a1]:
	
			dic[a1][a2]=0

f.close()

for a1 in sorted(dic.keys()):

	for a2 in sorted(dic[a1].keys()):

		print a1 + '\t' + a2 

