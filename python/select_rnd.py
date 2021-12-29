#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import regex as re


input_fasta = sys.argv[1]
rnd = sys.argv[2]
output_fasta = sys.argv[3]
my_records = []


for seq_record in SeqIO.parse(input_fasta, "fasta"):
    try:
        re.search(f'.*{rnd}.*', seq_record.id).group()
        rec = SeqRecord(seq_record.seq, id=seq_record.id, description='')
        my_records.append(rec)
    except:
        pass

SeqIO.write(my_records, output_fasta, 'fasta')