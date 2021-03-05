#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys



input_fasta = sys.argv[1]
gene_file = sys.argv[2]
output_fasta = sys.argv[3]
gene_lst = []
my_records = []

with open(gene_file) as genes:
    genes = genes.readlines()
    for gene in genes[1:]:
        gene = gene.strip().split(',')
        gene_lst.append(gene[1].strip('"'))

#print(gene_lst)

for seq_record in SeqIO.parse(input_fasta, "fasta"):
    if seq_record.id in gene_lst:
       # print(seq_record.id)
        rec = SeqRecord(seq_record.seq, id=seq_record.id, description='')
        my_records.append(rec)

SeqIO.write(my_records, output_fasta, 'fasta')