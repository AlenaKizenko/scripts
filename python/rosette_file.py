#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]

names_lst = [seq_record.id for seq_record in SeqIO.parse(input_file, 'fasta')]

with open(output_file, 'w') as file:
    for name in names_lst:
        fam_name = name.strip().split('_')[0]
        file.write(f'{name} {fam_name}\n')