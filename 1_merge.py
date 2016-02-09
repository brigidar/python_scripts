
#########################################################################################
#											#
# Name	      :	1_merge.py								#
# Version     : 1.4									#
# Project     : sort & merge SNP tables							#
# Description : Script to merge snp tables and recalculate syn, query_aa, query_codon.		#
# Author      : Brigida Rusconi								#
# Date        : April 29, 2015							#
#											#
#########################################################################################


#!/usr/bin/env python
# The merged table doesn't have dn/ds, snps/gene length, and transistion/transversion columns these will be calculated after the sorting script

#------------------------------------------------------------------------------------------

import argparse, os, sys, csv
import pandas
import pdb
import glob
import Bio
from pandas import *
from glob import glob
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
#------------------------------------------------------------------------------------------

#output file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="combined snp table and add csv extension to name")
parser.add_argument('-v', '--table', help="translation table selection", default=11)

args = parser.parse_args()
output_file = args.output
v=args.table
g=int(v)
#------------------------------------------------------------------------------------------
#calls up all txt files in the current directory

files = glob('*.txt')

# concatenates them all together and reads molecule and refpos as indexes
data = [read_csv(f, sep ='\t', header=0, dtype=unicode) for f in files]
data = [d.T for d in data]
d = concat(data, keys=files)
# transpose the table to have the correct headers
df = d.T
df=df.fillna('--')
#pdb.set_trace()

#pdb.set_trace()
#------------------------------------------------------------------------------------------


#redo the syn/nsyn/intergenic calling
count_qbase=df.columns.values
#condensate codon calling
codonindexes=[]
for i, v in enumerate(count_qbase):
    if 'query_codon' in v[1]:
        codonindexes.append(i)
codon=df.iloc[:,codonindexes]
pdb.set_trace()
# create list with all the unique values for the codons
CC=[]
posc=[]
for i in range(0,codon.index.size):
    posc=[n for n in unique(codon.iloc[i,:])[:]]
    CC.append(posc)

##join together all unique items with / as separator
CCS=[]

for i in CC:
    if len(i)==1:
        CCS.append(i[0])
    else:
        CCS.append('/'.join(set([n for n in '/'.join(i).split('/') if n != '--'])))

# second list where / are split up as list



#get qbase columns and ref_codon column
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v[1]:
        qindexes.append(i)
qbas=df.iloc[:,qindexes]
ref_cod=df.iloc[:,qindexes[-1]+7]
ref_codon=[]
for i in ref_cod:
    ref_codon.append(i)

#refbase column
df3=df.iloc[:,3].T
#refbase list
ref_base=[]
for i in df3:
    ref_base.append(i)
#pdb.set_trace()
#get position in codon
cod_pos=[]
for i,v in enumerate(ref_codon):
        if ref_codon[i] == CCS[i] :
            cod_pos.append(v)
        elif ref_codon[i][0]!=CCS[i][0]:
            cod_pos.append(0)
        elif ref_codon[i][1]!=CCS[i][1]:
            cod_pos.append(1)
        elif ref_codon[i][2]!= CCS[i][2]:
            cod_pos.append(2)

#get alleles for each snp
snp_nb=[]
for i in range (0,qbas.index.size):
    snp_uniq=[n for n in unique(qbas.iloc[i,:])[:] if n!= 'No Hit' and n!='indel']
    snp_nb.append(snp_uniq)
snp3=[]
for i,n in enumerate(snp_nb):
    snp4=[]
    for t in n:
        snp4+= t.split('/')
    snp3.append(snp4)


#check if codon is still needed for the snp present
fin_codon=[]
cod=[]
snp2=[]
for i,v in enumerate(CCS):
    cod=v.split('/')
    snp2=[n for n in snp3[i]]
    in_codon=[]
    if len(cod) == 1:
        if cod[0]==ref_codon[i]:
            fin_codon.append('--')
        else:
            fin_codon.append(cod[0])
    elif len(cod)>1:
        for c in cod:
            if c == ref_codon[i]:
                in_codon.append('--')
            elif c[cod_pos[i]] in snp2 or c=='indel':
                in_codon.append(c)
        fin_codon.append('/'.join(in_codon))




# with the set function it only keeps unique values because sometimes there is duplicates in the codons after merging
#------------------------------------------------------------------------------------------

bases=['A','C','G','T']
AA=[]

#translate each codon to the corresponding aa using bipython
for i in fin_codon:
    codons = [n for n in i.split('/') if n != '--' and n != 'indel']
    
# looks if there is only one empty one or an indel
    if len(codons) == 0:
        
        AA.append('--')
        
    else:
        
        aa = []
#makes an amino acid out of each codon the if all() statement checks that all nucleotides are ATGC. If there is anything else the codon will not be translated and will be identified as unknown.
        for p in codons:
            if all(z in bases for z in p)==True:
                aa.append(str(Seq(p, generic_dna).translate(table=g)))
            else:
                aa.append('Unknown')
        AA.append('/'.join(aa))


#------------------------------------------------------------------------------------------


# find columns that have qbase in it and add them together with the position and the refbase and molecule name




#adds refbase before qbase hits

df4=concat([df3,qbas], axis=1)


#------------------------------------------------------------------------------------------


# compare query amino acids to reference amino acid to get syn to nsyn and intergenic
refaa=[]
df2=df.reset_index()
for n,v in enumerate(df2.index):
    refaa.append(df2.iloc[n,len(qindexes)+11])

syn=[]
for i,v in enumerate(AA):
    amino = [n for n in AA[i].split('/') if n != '--' ]
    if len(amino) == 0:
        syn.append('intergenic')
    else:
        bb=[]
        for a in amino:
            if a==refaa[i]:
                bb.append('SYN')
            else:
                bb.append('NSYN')
        syn.append('/'.join(bb))

SYN=DataFrame(syn, index=df.index.values, columns=['syn?'])
df5=concat([SYN,df4], axis=1, join_axes=[df4.index])
df5_1=df.iloc[:,0:2]
df5_2=concat([df5_1,df5], axis=1, join_axes=[df5.index])

#------------------------------------------------------------------------------------------
# add transition and transversion information function developed by Mando Rodriguez

def get_trans_state(base, hit):
    
    if base == hit:
        return "--"
    elif base == 'A' and hit == 'G':
        return "transition"
    elif base == 'G' and hit == 'A':
        return "transition"
    elif base == 'C' and hit == 'T':
        return "transition"
    elif base == 'T' and hit == 'C':
        return "transition"
        
    elif base == 'A' and hit == 'C':
        return "transversion"
    elif base == 'C' and hit == 'A':
        return "transversion"
        
    elif base == 'A' and hit == 'T':
        return "transversion"
    elif base == 'T' and hit == 'A':
        return "transversion"
        
    elif base == 'C' and hit == 'G':
        return "transversion"
    elif base == 'G' and hit == 'C':
        return "transversion"
        
    elif base == 'G' and hit == 'T':
        return "transversion"
    elif base == 'T' and hit == 'G':
        return "transversion"
        
    else:
        return "Error [base:%s, hit:%s]" % (base, hit)


trs_trv=[]

for i,v in enumerate(ref_base):
    m=[n for n in snp3[i]]
    p=[]
    if len(m)>0:
        for t in m:
            p.append(get_trans_state(v,t))
        trs_trv.append('/'.join(p))
    else:
        trs_trv.append('--')

trs_trv2=[]
for i,v in enumerate(trs_trv):
    m=v.split('/')
    p=[]
    if len(m)>1:
        for t in m:
            if t !='--':
                p.append(t)
        trs_trv2.append('/'.join(p))
    else:
        trs_trv2.append(v)


# columns with refcodon refbase gene_name etc.
df6=df.iloc[:,(qindexes[-1]+1):(qindexes[-1]+9)]
#column with product


prod=df.iloc[:,df.columns.size-1]
prod2=[]
for n in prod:
    prod2.append(n)

df7=concat([df5_2,df6], axis=1, join_axes=[df5.index])


#------------------------------------------------------------------------------------------



#get num_hits columns and
nindexes=[]
for i, v in enumerate(count_qbase):
    if 'num_hits:' in v[1]:
        nindexes.append(i)
num=df.iloc[:,nindexes].T.reset_index(drop=True, level=0).T


#------------------------------------------------------------------------------------------

#column names have the text file name integrated in it and it needs to be removed
df8=df7.T.reset_index()

col_name=[]
for i in range(0,df8.index.size):
    col_name.append(df8.iloc[i,0][1])

col_name[2]='syn?'
COL=DataFrame(col_name, index=df8.index.values)
final=concat([COL,df8.iloc[:,1:]], axis=1)
final.set_index(final.iloc[:,0], inplace=True)
final=final.iloc[:,1:].T
final.insert(final.columns.size, 'query_codon', fin_codon)
final.insert(final.columns.size, 'query_aa', AA)
final2=concat([final,num],axis=1, join_axes=[final.index])
final.insert((final.columns.size), 'pos_in_codon', cod_pos)
final.insert(final.columns.size, 'transition/transversion', trs_trv2)
final.insert(final.columns.size, 'product', prod2)
final3=final.fillna('--')
#pdb.set_trace()
#------------------------------------------------------------------------------------------

#
final2=final2.astype('str').replace({'gene_name':{'intergenic':'None'}})
final2=final2.astype('str').replace({'query_aa':{'*':'Stop'}})
final2=final2.astype('str').replace({'ref_aa':{'*':'Stop'}})
final2.insert(final2.columns.size,'product', prod2)
with open('circleator_format.txt','w') as output:
    final2.to_csv(output, sep='\t', index=False)

# save new tab delimited table

with open(output_file,'w') as output:
 final3.to_csv(output, sep='\t', index=False)



