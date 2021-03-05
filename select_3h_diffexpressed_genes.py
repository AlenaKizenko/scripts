#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
import sys
import re


input_fasta = sys.argv[1]
gff_file = sys.argv[2]
gene_lst = []


for seq_record in SeqIO.parse(input_fasta, "fasta"):
    gene_lst.append(seq_record.id)


with open(gff_file) as genes:
    genes = genes.readlines()
    for line in genes:
        line = line.strip().split('\t')
        repid = re.search('(?<=Target=).+?(?= )', line[8]).group()
        if repid in gene_lst:
            print(
                f'{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\t{line[8]}\t')
