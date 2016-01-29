#########################################################################################
#											#
# Name	      :	coverage.py								#
# Version     : 1.0									#
# Project     : coverage of assembly							#
# Description : Script to verify coverage of assembly		#
# Author      : Brigida Rusconi								#
# Date        : May 18, 2015							#
#											#
#########################################################################################
import argparse, os, sys, csv, IPython
import pandas
import pdb
from pandas import *

parser = argparse.ArgumentParser()
parser.add_argument('-q', '--qc', help="mpileup file")
parser.add_argument('-c', '--coverage', help="coverage")

args = parser.parse_args()
cov=args.coverage
qc=args.qc

#reads in qc value file
df2 =read_csv(qc,sep='\t', header=None, names=['molecule','refpos','base','coverage','read','quality'])
df2=df2.set_index(['molecule','refpos'])
print " mpileup loaded"
print df2.index.size

df3=df2[df2['coverage']>=int(cov)]
print "removed low coverage"

print df3.index.size

