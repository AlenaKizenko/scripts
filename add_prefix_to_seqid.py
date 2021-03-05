#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# add prefix to sequence ids

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import statistics

def prefix_adding(filename, output_folder, prefix):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(filename)))[0]
    target_records = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        rec = SeqRecord(seq_record.seq, id = f'{prefix}_{seq_record.id}', description = '')
        target_records.append(rec)
    SeqIO.write(target_records, f'{output_folder}/{base_name}_prefix.fasta', 'fasta')

result = prefix_adding(sys.argv[1], sys.argv[2], sys.argv[3])