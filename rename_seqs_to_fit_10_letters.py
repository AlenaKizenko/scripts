#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# rename fasta ids to transform in Phylip
# arg1 = input filename, arg2 = output path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def split_seqid(filename, output_folder):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(filename)))[0]
    target_records = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        lst = seq_record.id.split('_')
        type_rep = lst[0]
        type_rnd = lst[2].split('-')[1]
        type_fam = lst[3].split('-')[1]
        type_num = lst[4]
        rec = SeqRecord(seq_record.seq, id=f'{type_rep[0]}_{type_rnd}_{type_fam}_{type_num}', description='')
        target_records.append(rec)
    SeqIO.write(target_records, f'{output_folder}/{base_name}_renamed.fasta', 'fasta')

result = split_seqid(sys.argv[1], sys.argv[2])