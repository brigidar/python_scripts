#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	1_merge_2.py								#
# Version     : 1.6
# Project     : sort & merge SNP tables							#
# Description : Script to merge snp tables and recalculate syn, query_aa, query_codon.		#
# Author      : Brigida Rusconi								#
# Date        : October 20, 2015							#
#											#
#########################################################################################



# The merged table doesn't have dn/ds, snps/gene length, and transistion/transversion columns these will be calculated after the sorting script

#------------------------------------------------------------------------------------------

import argparse, os, sys, csv
import pdb
import glob
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
#------------------------------------------------------------------------------------------


#redo the syn/nsyn/intergenic calling
count_qbase=df.columns.values
#condensate codon calling
codonindexes=[]
for i, v in enumerate(count_qbase):
    if 'query_codon' in v[1]:
        codonindexes.append(i)
codon=df.iloc[:,codonindexes]

# create list with all the unique values for the codons
CC=[]
posc=[]
for i in range(0,codon.index.size):
    posc=[n for n in unique(str(codon.iloc[i,0]).split('/'))[:]]
    CC.append(posc)
#pdb.set_trace()
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


for i, v in enumerate(count_qbase):
    if 'ref_codon' in v[1]:
        ref_codon=df.iloc[:,i]



#refbase list
for i, v in enumerate(count_qbase):
    if 'refbase' in v[1]:
        ref_base=df.iloc[:,i]

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

#df4=concat([df3,qbas], axis=1)


#------------------------------------------------------------------------------------------



# compare query amino acids to reference amino acid to get syn to nsyn and intergenic

for i, v in enumerate(count_qbase):
    if 'ref_aa' in v[1]:
        refaa=df.iloc[:,i]


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



SYN=DataFrame(syn, index=df.index.values, columns=['syn/nsyn/intergenic'])

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
#df6=df.iloc[:,(qindexes[-1]+1):(qindexes[-1]+9)]
#column with product


prod=df.iloc[:,df.columns.size-1]
prod2=[]
for n in prod:
    prod2.append(n)



#------------------------------------------------------------------------------------------



#get num_hits columns and
nindexes=[]
for i, v in enumerate(count_qbase):
    if 'num_hits:' in v[1]:
        nindexes.append(i)
num=df.iloc[:,nindexes].T.reset_index(drop=True, level=0).T


#------------------------------------------------------------------------------------------

#column names have the text file name integrated in it and it needs to be removed
df2=df.T.reset_index(drop=True, level=0).T


final=df2.iloc[:,0:2]

final.insert(final.columns.size, 'refbase', ref_base)
final2=concat([final,df2.iloc[:,qindexes[0]:(qindexes[0]+9)]], axis=1, join_axes=[final.index])
final2.insert(final2.columns.size, 'query_codon', fin_codon)
final2.insert(final2.columns.size, 'query_aa', AA)

# format for circelator
final3=concat([final2,num],axis=1, join_axes=[final.index])
final3.insert(final3.columns.size, 'product', prod2)
final3.insert(2, 'syn?', syn)
final3.fillna('--')

# format for merged table
final2.insert(2, 'syn/nsyn/intergenic', syn)
final2.insert(final2.columns.size, 'transition/transversion', trs_trv2)
final2.insert((final2.columns.size), 'pos_in_codon', cod_pos)
final2.insert(final2.columns.size, 'product', prod2)





#pdb.set_trace()
#------------------------------------------------------------------------------------------

#
final3=final3.astype('str').replace({'gene_name':{'intergenic':'None'}})
final3=final3.astype('str').replace({'query_aa':{'*':'Stop'}})
final3=final3.astype('str').replace({'ref_aa':{'*':'Stop'}})
with open('circleator_format.txt','w') as output:
    final3.to_csv(output, sep='\t', index=False)

# save new tab delimited table

with open(output_file,'w') as output:
 final2.to_csv(output, sep='\t', index=False)



