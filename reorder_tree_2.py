#########################################################################################
#											#
# Name	      :	reorder_tree_2.py							#
# Version     : 1.0									#
# Project     : SNP Genotyper							#
# Description : Script to reorder snp table format by phylogeny.
# Author      : Brigida Rusconi								#
# Date        : June 10, 2015							#
#											#
#########################################################################################
# for isin information
#http://pandas.pydata.org/pandas-docs/stable/indexing.html

#!/usr/bin/env python


import argparse, os, sys, csv, IPython
import pandas
import pdb
from pandas import *
from IPython import get_ipython
import matplotlib.pyplot as plt
from pandas.util.testing import assert_frame_equal
import Bio
from Bio import Phylo


#output and input file name to give with the script and arguments for cutoff
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="SNPs reordered according to tree")
parser.add_argument('-s', '--snp_table', help="circleator table")
parser.add_argument('-t', '--tree', help="newick formatted tree")



args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
tree1 = args.tree

#read in files as dataframe
df =read_csv(input_file, sep='\t', header=0)

#gets order of tree
def parse_tree(tree):
    names = []
    mytree = Phylo.read(tree, 'newick')
    tree_names = [ ]
    for clade in mytree.find_clades():
        if clade.name:
            names.append(clade.name)
    return names
#uses function parse_tree to get list of names
trees=parse_tree(tree1)

#get header with qbase
count_qbase=df.columns.values
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(v)
# reorder snp file based on tree
df2= DataFrame(index=df.index.values)
for name in trees:
    for i in qindexes:
        if i.split(':')[1]==name:
            df2=concat([df2,df.loc[:,i]], axis=1)

##get header with num_hits
#nindexes=[]
#for i, v in enumerate(count_qbase):
#    if 'num_hits:' in v:
#        nindexes.append(v)
#
#df3= DataFrame(index=df.index.values)
#for name in trees:
#    for i in nindexes:
#        if i.split(':')[1]==name:
#            df3=concat([df3,df.loc[:,i]], axis=1)
#
#put back together
df4=concat([df.iloc[:,:4], df2], axis=1)

df5=df.iloc[:, len(qindexes)+4:]

df6=concat([df4,df5], axis=1)

#save reordered file
with open(output_file,'w') as output:
    df6.to_csv(output, sep='\t', index=False)






                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 


