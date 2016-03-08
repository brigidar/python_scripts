import argparse, os, sys
import pdb
from Bio import SeqIO


#argument to give to the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="filtered proteins")
parser.add_argument('-i', '--input', help="total proteins fasta")
parser.add_argument('-c', '--id', help="id list")


args = parser.parse_args()
output_file = args.output
input_file = args.input
id_num = args.id
id_list=[]
with open(id_num, 'r') as item:
     id_list=item.read().split('\n')




#open fasta file and select everything above cutoff:

seq_iterator = []
for record in SeqIO.parse(open(input_file, 'rU'), "fasta"):
    if ">%s" %record.id in id_list:
        seq_iterator.append(record)

# save in output
with open(output_file, 'w') as output_handle:
    SeqIO.write(seq_iterator, output_handle, "fasta")



