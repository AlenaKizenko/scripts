#!/usr/bin/env python3


import sys
import re


gff_file = sys.argv[1]
gene_lst = []


with open(gff_file) as genes:
    genes = genes.readlines()
    for line in genes:
        line = line.strip().split('\t')
        target = re.search('(?<=Target=).+?(?= )', line[8]).group()
        id = re.search('(?<=ID=).+?(?=;)', line[8]).group()
        col = line[8].strip().split(' ')
        #print(" ".join(col[1:4]))
        #print(col[1:4])
        #print(f'ID={target};Target={id} {" ".join(col[1:3])}')
        print(
            f'{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t'
            f'{line[5]}\t{line[6]}\t{line[7]}\tID={target};Target={id} {" ".join(col[1:3])}\n')