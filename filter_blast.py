#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	filter_blastn.py								#
# Version     : 0.1									#
# Project     : sort & merge SNP tables							#
# Description : find blsatn hits of IS elements in reference genome	#
# Author      : Brigida Rusconi								#
# Date        : March 8th, 2016							#
#											#
#########################################################################################
#for replacement of a given value with NaN
#http://stackoverflow.com/questions/18172851/deleting-dataframe-row-in-pandas-based-on-column-value

# to remove any symbol and replace it with nan
#http://stackoverflow.com/questions/875968/how-to-remove-symbols-from-a-string-with-python

# for isin information
#http://pandas.pydata.org/pandas-docs/stable/indexing.html

# for selecting rows that have an indel:
#http://stackoverflow.com/questions/14247586/python-pandas-how-to-select-rows-with-one-or-more-nulls-from-a-dataframe-without



#------------------------------------------------------------------------------------------


import argparse, os, sys, csv
import pdb
import numpy as np
from pandas import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
#------------------------------------------------------------------------------------------


#output and input file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="sorted filtered table for fasta")
parser.add_argument('-b', '--blastn', help="blastn input")


args = parser.parse_args()
output_file = args.output
input_file = args.blastn
#------------------------------------------------------------------------------------------


#read in file as dataframe
df =read_csv(input_file,sep='\t', header=None, dtype=object)
#need to fill na otherwise it cannot do boolean operations

#pdb.set_trace()
for i in range(0,df.index.size):
    if float(df.iloc[i,3])/float(df.iloc[i,23])<0.5:
        df.iloc[i,3]=0

df2=df[df.iloc[:,3]>0]
#pdb.set_trace()
df3=concat([df2.iloc[:,0],df2.iloc[:,6:8]], axis=1)
new_index=[]
new_index2=[]
for i in range(0,df3.index.size):
    if df3.iloc[i,1]>df3.iloc[i,2]:
        new_index.append(df3.index[i])
    else:
        new_index2.append(df3.index[i])


#pdb.set_trace()
#inverted matches (start bigger than stop)
if len(new_index)>0:
    df4=df3[df3.index.isin(new_index)]
#pdb.set_trace()
    cols = df4.columns.tolist()
#pdb.set_trace()
    cols=[cols[0]]+cols[-1:]+[cols[1]]
    df4=df4[cols]
#normal order
    df5=df3[df3.index.isin(new_index2)]

    df6=concat([df5,df4],axis=1)
else:
    df6=df5
pdb.set_trace()
#------------------------------------------------------------------------------------------




#save total file for plotting -t option
with open(output2_file,'w') as output2:
    df6.to_csv(output2, sep='\t', index=False, header=None)


