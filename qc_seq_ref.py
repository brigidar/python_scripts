
#########################################################################################
#											#
# Name	      :	qc_seq_ref.py								#
# Version     : 1.0									#
# Project     : qc of verified snps							#
# Description : Script to verify coverage of refbase		#
# Author      : Brigida Rusconi								#
# Date        : April 29, 2015							#
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

args = parser.parse_args()
input_file = args.snp_table
output=args.output
cov=args.coverage
qc=args.qc

#reads in snp table
df =read_csv(input_file,sep='\t')
df=df.set_index(['molecule','refpos'])
print " snp table loaded"
#reads in qc value file
df2 =read_csv(qc,sep='\t', header=None, names=['molecule','refpos','base','coverage','quality','6'])
df2=df2.set_index(['molecule','refpos'])
print " mpileup loaded"

#reindex qc file to snp pos
df3=df2[df2.index.isin(df.index)]
print " reindexed"

#keep only those with a value above coverage given
df4=df3[df3['coverage']>=int(cov)]
print "removed low coverage"
#reindex according to coverage of query genome
df5=df[df.index.isin(df4.index)]
print "reindex according to coverage"
with open(output,'w') as output:
    df5.to_csv(output, sep='\t')
