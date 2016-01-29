
#########################################################################################
#											#
# Name	      :	qc_seq.py								#
# Version     : 1.0									#
# Project     : qc of verified snps							#
# Description : Script to verify coverage of hit bases		#
# Author      : Brigida Rusconi								#
# Date        : April 30, 2015							#
#											#
#########################################################################################
import argparse, os, sys, csv, IPython
import pandas
import pdb
from pandas import *

parser = argparse.ArgumentParser()
parser.add_argument('-q', '--qc', help="mpileup file")
parser.add_argument('-s', '--snp_table', help="snp table for snp distribution")
parser.add_argument('-o', '--output', help="output")
parser.add_argument('-c', '--coverage', help="coverage")
parser.add_argument('-d', '--query_qc', help="positions in blast for hits")
parser.add_argument('-b', '--column1', help="column1", type=int)
parser.add_argument('-g', '--column2', help="column2", type=int)

args = parser.parse_args()
input_file = args.snp_table
output=args.output
cov=args.coverage
query_qc=args.query_qc
b=args.column1
g=args.column2
qc=args.qc
#reads in snp table
df =read_csv(input_file,sep='\t')
df=df.set_index(['molecule','refpos'])
print " snp table loaded"
#reads in qc value file
df2 =read_csv(qc,sep='\t', header=None, names=['molecule','refpos','base','coverage','read','quality'])
df2=df2.set_index(['molecule','refpos'])
print " mpileup loaded"

#reads in file with query qc
df3=read_csv(query_qc, sep='\t')
print "blast fetch loaded"
headers=list(df3.columns.values)


df4=df3.set_index(headers[b:g])
#reindex qc file to snp pos
df5=df2[df2.index.isin(df4.index)]
print " reindexed"
#keep only those with a value above coverage given
df6=df5[df5['coverage']>=int(cov)]
print "removed low coverage"
#reindex according to coverage of query genome
df4=df4[df4.index.isin(df6.index)]
df4.reset_index()
new=df4.set_index(['molecule','refpos'])
print "reindex according to coverage"
#reindex to high coverage positions
df7=df[df.index.isin(new.index)]

with open(output,'w') as output:
    df7.to_csv(output, sep='\t')
