#########################################################################################
#											#
# Name	      :	circleator-format.py								#
# Version     : 1.0								#
# Project     : SNP table for circleator 							#
# Description : Script tohave the right format		#
# Author      : Brigida Rusconi								#
# Date        : March 30, 2015							#
#											#
#########################################################################################


#!/usr/bin/env python
# if the NC number has a .1 at the end just do a sed to remove it after the format changes

#------------------------------------------------------------------------------------------

import argparse, os, sys, csv, IPython
import pandas
import pdb
import glob
import Bio
from pandas import *
from IPython import get_ipython
from glob import glob
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
#------------------------------------------------------------------------------------------

#output file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="tabel for circleator")
parser.add_argument('-s', '--input', help="filtered table", default=11)

args = parser.parse_args()
output_file = args.output
input=args.input

df =read_csv(input, sep='\t', index_col=None, dtype=unicode)


# table for circleator
filler=[0 for i in range(0,df.index.size)]
filler2=['nucmer'for i in range(0,df.index.size)]
#pdb.set_trace()
df.insert((df.columns.size-2), 'properties', filler2)
df.insert((df.columns.size-2), 'Num_No_Hit', filler)
df.insert((df.columns.size-2), 'Homoplasy', filler)
df=df.astype('str').replace({'syn?':{'intergenic':'NA'}})
#pdb.set_trace()
first=df.loc[:,'molecule':'gene_name']
second=df.loc[:, 'gene_start':'Homoplasy']
df2=concat([first,df['product']], axis =1)
df3=concat([df2,second],axis=1)
#pdb.set_trace()
df4=df3.T.reset_index()
col_name=[]
for i in range(0,df4.index.size):
    col_name.append(str(df4.iloc[i,0]))
for i,v in enumerate(col_name):
    if v=='gene_end':
        col_name[i]='gene_stop'
col_name2=[]
for i,v in enumerate(col_name):
    col=v.split('.')
    if len(col)==1:
        col_name2.append(v)
    else:
        col_name2.append(col[0])

col_name3=[]
for i,v in enumerate(col_name2):
    col=v.split('qbase:')
    if len(col)==1:
        col_name3.append(v)
    else:
        col_name3.append(col[1])
COL=DataFrame(col_name3, index=df4.index.values)

df5=concat([COL,df4.iloc[:,1:]], axis=1)
df5.set_index(df5.iloc[:,0], inplace=True)
df5=df5.iloc[:,1:].T
#pdb.set_trace()


with open(output_file,'w') as output:
    df5.to_csv(output, sep='\t', index=False)
