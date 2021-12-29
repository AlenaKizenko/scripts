#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


input_fasta = sys.argv[1]

for seq_record in SeqIO.parse(input_fasta, "fasta"):
    rec = SeqRecord(seq_record.seq, id=seq_record.id, description='')  # create SeqRecord object
    SeqIO.write(rec, f'{seq_record.id}.fasta', 'fasta')