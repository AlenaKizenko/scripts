#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


paf_file = sys.argv[1]
input_fasta = sys.argv[2]
align_out = sys.argv[3]


with open(paf_file) as paf:
    paf = paf.readlines()
    paf_dict = {}
    cnt = 1
    for line in paf:
        line = line.strip().split('\t')
        paf_dict[cnt] = []
        paf_dict[cnt].append(line[2])
        paf_dict[cnt].append(line[3])
        cnt += 1
        

rec_lst = []
seq_record = SeqIO.read(input_fasta, "fasta")

for i in range(1, 10):
    coord1 = int(paf_dict[i][0])
    coord2 = int(paf_dict[i][1])
    seq_record_seq_before = seq_record.seq[:coord1]
    rec_before = SeqRecord(seq_record_seq_before, id=f'Alignment 0:{coord1}', description='')
    rec_lst.append(rec_before)
    #seq_record_seq_after = seq_record.seq[coord2:]
    #rec_after = SeqRecord(seq_record_seq_after, id=f'Alignment {coord2}:end', description='')
    #rec_lst.append(rec_after)

SeqIO.write(rec_lst, f'{align_out}', 'fasta')