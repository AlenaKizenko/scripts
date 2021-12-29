#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
from Bio.Seq import Seq

file = sys.argv[1]
output_folder = sys.argv[2]
exp = []
rst = []
target_records = []

for seq_record in SeqIO.parse(file, "fasta"):
    if seq_record.id[0:3] == 'exp':
        exp.append(str(seq_record.seq))
    elif seq_record.id[0:3] == 'rst':
        rst.append(str(seq_record.seq))



exp_str = ''.join(exp)
rec_exp = SeqRecord(Seq(exp_str), id='LINE_CR1_expansion', description='')
target_records.append(rec_exp)


rst_str = ''.join(rst)
rec_rst = SeqRecord(Seq(rst_str), id='LINE_CR1_rest', description='')
target_records.append(rec_rst)

base_name = os.path.splitext(os.path.basename(os.path.normpath(file)))[0]
SeqIO.write(target_records, f'{output_folder}/{base_name}_conc.fasta', 'fasta')

