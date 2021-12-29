#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def split_seqid(new, old, rt1_file, output_folder):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(rt1_file)))[0]
    target_records = []
    new_lst = [line.strip() for line in open(new)]
    old_lst = [line.strip() for line in open(old)]
    for seq_record3 in SeqIO.parse(rt1_file, "fasta"):
        seq_record3.id = seq_record3.id.split('.p')[0]
        if seq_record3.id in new_lst:
            rec = SeqRecord(seq_record3.seq, id=f'new_{seq_record3.id}', description='')
            target_records.append(rec)
        elif seq_record3.id in old_lst:
            rec = SeqRecord(seq_record3.seq, id=f'old_{seq_record3.id}', description='')
            target_records.append(rec)
        else:
            pass
    SeqIO.write(target_records, f'{output_folder}/{base_name}_prefix.fasta', 'fasta')


result = split_seqid(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])