#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def split_seqid(filename, output_folder):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(filename)))[0]
    cnt = 0
    target_records = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        cnt += 1
        seq_record.id = seq_record.id.split(' ')[0]
        seq_record.id = seq_record.id.replace("(-)", "_minus")
        seq_record.id = seq_record.id.replace("(+)", "_plus")
        rec = SeqRecord(seq_record.seq, id = seq_record.id, description = '')
        target_records.append(rec)
    SeqIO.write(target_records, f'{output_folder}/{base_name}_renamed.fasta', 'fasta')
    
result = split_seqid(sys.argv[1], sys.argv[2])