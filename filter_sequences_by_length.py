#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# filter input sequences by sequence length

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import statistics

def mean_filtering(filename, output_folder):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(filename)))[0]
    len_lst = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]
    mean_len = statistics.mean(len_lst)
    target_records = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        if len(seq_record) >= mean_len:
            rec = SeqRecord(seq_record.seq, id = seq_record.id, description = '')
            target_records.append(rec)
    SeqIO.write(target_records, f'{output_folder}/{base_name}_filtered.fasta', 'fasta')

result = mean_filtering(sys.argv[1], sys.argv[2])