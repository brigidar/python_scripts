
#########################################################################################
#											#
# Name	      :	5_sorting_genotyper2_stop.py								#
# Version     : 1.1									#
# Project     : SNP table downstream analysis						#
# Description : Script to summarize the genotyper output and classify according to groups etc.		#
# Author      : Brigida Rusconi								#
# Date        : April 1, 2015							#
#											#
#########################################################################################

import argparse, os, sys, csv
import pandas
import pdb
from pandas import *
from pandas.util.testing import assert_frame_equal

#------------------------------------------------------------------------------------------------------------

#output and input file name
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="genotype group summary")
parser.add_argument('-s', '--genotyper_table', help="genotyper table to summarize", default="snp_genotype_out.txt")

parser.add_argument('-r', '--refgenome', help="refgenome name")

args = parser.parse_args()
output_file = args.output
input_file = args.genotyper_table
ref_genome= args.refgenome

#------------------------------------------------------------------------------------------------------------


#read in files as dataframe
df =read_csv(input_file, sep='\t', header=0, dtype=str)


#------------------------------------------------------------------------------------------------------------

##get sums of total patterns, SYN, NSYN, Intergenic according to PI, NI and multiallelic state

total=[df.index.size]


# get total of PI NI and add to total snps
pl=[]
nl=[]
for name, group in df.groupby('Informative'):
    nl.append(name)

for i in df.groupby('Informative').size():
    pl.append(i)
#------------------------------------------------------------------------------------------------------------


#total amount of stop
stop=[str(df['query_aa'].str.contains('Stop').sum())]
#------------------------------------------------------------------------------------------------------------

#total amount of hypothetical
hypo=[str(df['product'].str.contains('hypothetical').sum())]
#------------------------------------------------------------------------------------------------------------

#total amount of multiallelic positions
multi=[str(df['query_codon'].str.contains('/').sum())]
#------------------------------------------------------------------------------------------------------------

# counting number of transitions transversions:

ti_n=[]
ti_g=[]
for name, group in df.groupby('transition/transversion'):
    ti_n.append(name)
for i in df.groupby('transition/transversion').size():
    ti_g.append(i)

# make list with transition and transversion only to add to it the others with multiallelic state

ti_n2=[]
ti_g2=[]
for i,v in enumerate(ti_n):
    m=v.split('/')
    if len(m)==1:
        ti_n2.append(v)
        ti_g2.append(ti_g[i])



for i,v in enumerate(ti_n):
    m=v.split('/')
    if len(m)>1:
        for t in m:
            if t==ti_n2[0]:
                ti_g2[0]=ti_g2[0] + ti_g[i]
            elif t==ti_n2[1]:
                ti_g2[1]=ti_g2[1] + ti_g[i]



#------------------------------------------------------------------------------------------------------------
#first row total values
tot=total+pl+stop+hypo+ti_g2+multi
#------------------------------------------------------------------------------------------------------------

#counting number of syn/nsyn
other_syn=DataFrame(df.groupby('syn?').size())


#------------------------------------------------------------------------------------------------------------

#distribution of informative non-informative
pi_dis= DataFrame(df.groupby(['syn?', 'Informative']).size()).reset_index(level=1)
test=DataFrame(index=pi_dis.index)
for i,group in pi_dis.groupby('Informative'):
    test=concat([test,group[0]], axis=1)
test2=concat([other_syn,test], axis=1)

#------------------------------------------------------------------------------------------------------------

# distribution of stops
si_dis= DataFrame(df.groupby(['query_aa','syn?']).size()).reset_index(level=0)
tes=si_dis[si_dis['query_aa'].str.contains('Stop')]
si_dis2=tes.reset_index().groupby('syn?').sum()
test3=concat([test2,si_dis2], axis=1, join_axes=[test.index])
#------------------------------------------------------------------------------------------------------------


# distribution of hypothetical
hy_dis= DataFrame(df.groupby(['product', 'syn?']).size()).reset_index(level=0)
hy2=hy_dis[hy_dis['product'].str.contains('hypothetical')]
hy3=hy2.reset_index().groupby('syn?').sum()
test4=concat([test3,hy3], axis=1, join_axes=[test.index])
#------------------------------------------------------------------------------------------------------------


#distribution of transition transversion
ti_dis=DataFrame(df.groupby(['syn?', 'transition/transversion']).size()).reset_index()
ti=DataFrame(index=ti_dis.index)

tl=[]
r=[]
rl=[]

#get one value for transition and one for transversion

for i,group in ti_dis.groupby('syn?'):
    tl.append([n for n in group['transition/transversion'].values])
    rl.append([n for n in group.iloc[:,2].values])#[n for n in
trans=[0 for i in range(0,len(rl))]
trasv=[0 for i in range(0,len(rl))]
for i,m in enumerate(tl):
    for t,c in enumerate(m):
        p=c.split('/')
        for d in p:
            if d=='transition':
                trans[i]+= rl[i][t]
        
            elif d=='transversion':
                trasv[i]+= rl[i][t]



test4.insert(test4.columns.size, 'transition', trans)
test4.insert(test4.columns.size, 'transversion', trasv)





#------------------------------------------------------------------------------------------------------------


# distribution of multiallelic
mu_dis= DataFrame(df.groupby(['query_codon', 'syn?']).size()).reset_index(level=0)
mu2=mu_dis[mu_dis['query_codon'].str.contains('/')]
mu3=mu2.reset_index().groupby('syn?').sum()


#get ditsribution of the multiallelic states in the intergenic regions
mu_int=df.groupby('syn?').get_group('intergenic')
counter=0
for i,n in enumerate(mu_int['snp_total']):
    if len(n)==2:
        counter+=1
    elif len(n)==3:
        counter+=2
mu4=mu3.T
mu4.insert(mu4.columns.size, 'intergenic', counter)
mu5=mu4.T

#------------------------------------------------------------------------------------------------------------

test5=concat([test4,mu5], axis=1, join_axes=[test.index])
test5=test5.T
test5.insert(0,'total',tot)
col=['SNPs']+[n for n in nl]+['stop']+['hypothetical']+[n for n in ti_n2]+['multiallelic']
test5.insert(0, 'header',col)
test6 = test5.fillna('0').set_index('header').astype('float64')
# check back to string
#------------------------------------------------------------------------------------------------------------

#calculate percentage
perc=test6.copy(deep=True)
perc2=DataFrame(index=[test6.index])
#adds percentage column after each column first row is percentage compared to total other rows are percentage compared to subcategory (SYN, NSYN. etc) have to add number of multiallelic to have the right amount
for i in range(0,test6.columns.size):
    mol=[n for n in (test6.iloc[:,i] / test6.iloc[0,i]  * 100)]
    # since transition and transversion are separate for the multiallelic states the total has to be double for two state does not account for 3 state yet.
    mol1=[n for n in (test6.iloc[-3:-1,i] / (test6.iloc[0,i] + test6.iloc[-1,i])  * 100)]
    mol2=[test6.iloc[0,i]/test6.loc[test6.index[0],'total']*100]
    mol1.append(mol[-1])
    mol3=mol2+mol[1:-3]+mol1
#    pdb.set_trace()
    perc.insert(i+i+1,('% '+ str(test6.columns[i])),mol3)
    perc2.insert(i,('% '+ str(test6.columns[i])),mol3)
#------------------------------------------------------------------------------------------------------------

#save summary table with percentage
with open(output_file ,'w') as output:
    perc.to_csv(output, sep='\t')
#------------------------------------------------------------------------------------------------------------

#save summary only percentage with sorted index
with open("percentage_summary.txt" ,'w') as output:
    perc2.T.sort_index().to_csv(output, sep='\t')
#------------------------------------------------------------------------------------------------------------

#summary of all groups and positions and genome names #http://wesmckinney.com/blog/filtering-out-duplicate-dataframe-rows/

grouped = df.groupby('Group')
index = [gp_keys[0] for gp_keys in grouped.groups.values()]
unique = df.reindex(index).iloc[:,3:12].set_index('Group')
#groups the amount of syn nsyn intergenic for each one of the groups
ls1=df.groupby(['Group', 'syn?']).size()
#keeps only the group as an index
ls2=ls1.reset_index(level=1)
#change the name of the counting column
ls2.rename(columns={0:'count_snp'}, inplace=True)
#groups the amount of PI?NI for each one of the groups after dividing by SYN/NSYN
lsp=df.groupby(['Group', 'syn?','Informative']).size()
#keeps only the group as an index
lsp2=lsp.reset_index(level=[1,2])
#change the name of the counting column
lsp2.rename(columns={0:'count_PI'}, inplace=True)
#------------------------------------------------------------------------------------------------------------

#take the refpos that match that group and put them with the molecule name
refpos_list=[]
mol_list=[]
grouped=df.groupby(['Group', 'syn?',])
for item in grouped['refpos']:
    rf=[str(n) for n in item[1]]
    refpos_list.append('/'.join(rf))
for i in grouped['molecule']:
    mf=[str(n) for n in i[1]]
    #only unique value no duplicates
    mol_list.append(''.join(set(mf)))
ps=[]

for v,i in enumerate(mol_list):
         ps.append("%s:%s" % (i,refpos_list[v]) )

#concatenates based on the longer index with the repeated groups with the join_axes command
lst=concat([ls2, unique],axis=1, join_axes=[ls2.index])
ni_ref=[]
#------------------------------------------------------------------------------------------------------------


#add column for group that is specific to the reference genome. if the group has only one kind of mutation syn/intergenic, etc then it will call a string instead of a series so we need to check first the class of the loc calling. takes argument with ref genome name.
for i in lst.index:
    m=[]
    pat = ""
    
    #if lst.loc[i, 'Pattern'].__class__ == 'str':
    if isinstance(lst.loc[i, 'Pattern'], str):
    
        pat = lst.loc[i, 'Pattern']
    
    else:
    
        pat = lst.loc[i, 'Pattern'][0]
    
    
    for n in pat:
        m.append(int(n))

    if all(x==1 for x in m[1:]):
        ni_ref.append(ref_genome)
    else:
        ni_ref.append('--')
            
#------------------------------------------------------------------------------------------------------------


#remove old PI column
lst.insert(5, 'refpos split up',ps)
lst.insert(lst.columns.size,'ref_genome', ni_ref)
#------------------------------------------------------------------------------------------------------------


#save summary for group
with open("summary_groups.txt",'w') as output:
    lst.to_csv(output, sep='\t')












