#!/usr/bin/env python
##!/usr/local/anaconda/bin/python

#########################################################################################
#											#
# Name	      :	vcf_snp.py								#
# Version     : 0.3									#
# Project     : extract snp from vcf						#
# Description : Script to exctract snps		#
# Author      : Brigida Rusconi								#
# Date        : November 06 2015							#
#											#
#########################################################################################

# If the position is identical to the reference it does not print the nucleotide. I have to retrieve it from the ref column.


#------------------------------------------------------------------------------------------


import argparse, os, sys, csv, IPython
import pandas
import pdb
from pandas import *

#------------------------------------------------------------------------------------------


#output and input file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="snp tab")
parser.add_argument('-s', '--snp_table', help="vcf")
parser.add_argument('-p', '--positions', help="positions")



args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
input_file2=args.positions
#------------------------------------------------------------------------------------------


#read in file as dataframe
df =read_csv(input_file,sep='\t', dtype=object)
df3=read_csv(input_file2,sep='\t', dtype=object, names=['#CHROM','POS'])
#pdb.set_trace()
#get the positions that are missing
if df.index.size!= df3.index.size:
    df4=df3[~df3.POS.isin(df.POS)]
#pdb.set_trace()
    df4=df4.set_index('POS')
    df=df.set_index('POS')
    df2=concat([df,df4], axis=1)
#if there is no information for a loci the position is dropped so it fills the nan with N
    df5=df2.iloc[:, 2:4].fillna('N')
    df5=df5.reindex(df3.POS)
else:
    df5=df.iloc[:,3:5]

#------------------------------------------------------------------------------------------


#pdb.set_trace()

#------------------------------------------------------------------------------------------
#replaces the . with the nucleotide call of the reference also deals with multiallelic states calling them N
ref_list=[]
for i in range(0,df5.index.size):
    if df5.iloc[i,1]==".":
        ref_list.append(df5.iloc[i,0][0])
    elif len(df5.iloc[i,1])>1:
        ref_list.append('N')
    else:
        ref_list.append(df5.iloc[i,1][0])
#pdb.set_trace()
#
##------------------------------------------------------------------------------------------
#
#save file with output name for fasta -o option and removes header and index
with open(output_file,'w') as output:
    output.write('ALT' + '\t' + ''.join([str(i) for v,i in enumerate(ref_list)]))
##------------------------------------------------------------------------------------------