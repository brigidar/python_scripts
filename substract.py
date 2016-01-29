#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	substract.py								#
# Version     : 0.1								#
# Project     : normalize positions							#
# Description : normailze positions from snp to bed file		#
# Author      : Brigida Rusconi								#
# Date        : October 21, 2015							#
#											#
#########################################################################################


import argparse, os, sys, csv, IPython
import pandas
import pdb
from pandas import *



#output file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-s', '--input', help="table to transform")
parser.add_argument('-v', '--number', help="value to substract", default=1)
parser.add_argument('-o', '--output', help="output file", default="output.bed")

args = parser.parse_args()
input_file = args.input
v=args.number
output=args.output
g=int(v)

infile=open(input_file, 'r')
table=infile.readlines()
infile.close()

# skips first column with chr or locus_tag

table2=[]
for i,v in enumerate(table):
    line=[]
    for c,n in enumerate(v.split('\t')):
        if c==0:
            line.append(n)
        else:
            line.append(int(n)-g)
    table2.append(line)



with open(output,'w') as of:
    for i,v in enumerate(table2):
            of.write('\t'.join([str(i) for i in v])+"\n")