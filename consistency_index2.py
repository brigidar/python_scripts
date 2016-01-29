# filters the SNP according to the conisistency data that was obtained by paup. adds the conistency index to the snp table and gives two files according to a cutoff value for the CI.

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


#output and input file name to give with the script and arguments for cutoff
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="SNPs filtered according to CI")
parser.add_argument('-s', '--snp_table', help="snp table to filter out homoplasy")
parser.add_argument('-k', '--cutoff', help="CI cutoff value, default value is 0.6", type=int, default= 0.6)
parser.add_argument('-d', '--deleted', help = " Removed positions", default= "removed.txt")


args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
cutoff_val = args.cutoff
delt = args.deleted
#read in files as dataframe
df =read_csv(input_file, sep='\t', header=0)



# remove all lines that are below cutoff value
ci2 = df[df.CI.astype('float64') >= cutoff_val]

ci3 = df[df.CI.astype('float64') < cutoff_val]


#save ci above cutoff
with open(output_file,'w') as output:
    ci2.to_csv(output, sep='\t', index=False)

#save ci below cutoff
with open(delt,'w') as output:
    ci3.to_csv(output, sep='\t', index=False)





                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 


