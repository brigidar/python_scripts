
#########################################################################################
#											#
# Name	      :	reindex.py								#
# Version     : 1									#
# Project     : confirm snps							#
# Description : Script to reindex a bigger table according to snps selected in a smaller subset	#
# Author      : Brigida Rusconi								#
# Date        : June 10, 2015							#
#											#
#########################################################################################

#!/usr/bin/env python
#------------------------------------------------------------------------------------------


import argparse, os, sys, csv, IPython
import pandas
import pdb
import numpy as np
from pandas import *
from IPython import get_ipython
import matplotlib.pyplot as plt
from pandas.util.testing import assert_frame_equal
#------------------------------------------------------------------------------------------


#output and input file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="reduced file according to smaller index")
parser.add_argument('-s', '--snp_table', help="snp table to reindex")
parser.add_argument('-r', '--reindex', help="table with smaller index")


args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
reindex = args.reindex
#------------------------------------------------------------------------------------------


#read in file as dataframe
df =read_csv(input_file,sep='\t', dtype=object)
df=df.set_index(['molecule','refpos'])

df2 =read_csv(reindex,sep='\t', dtype=object)
df2=df2.set_index(['molecule','refpos'])


df3=df[df.index.isin(df2.index)]

with open(output_file,'w') as output:
    df3.to_csv(output, sep='\t')
