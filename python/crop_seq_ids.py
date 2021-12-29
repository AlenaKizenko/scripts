#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def split_seqid(rt1_file, output_folder):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(rt1_file)))[0]
    target_records = []
    for seq_record in SeqIO.parse(rt1_file, "fasta"):
        seq_record.id = seq_record.id.split('.p')[0]
        #fam = seq_record.id.split('#')[0]
        #clas = seq_record.id.split('#')[1]
        #rec = SeqRecord(seq_record.seq, id=f'{clas}_{fam}', description='')
        rec = SeqRecord(seq_record.seq, id=seq_record.id, description='')
        target_records.append(rec)
    SeqIO.write(target_records, f'{output_folder}/{base_name}_final.fasta', 'fasta')


result = split_seqid(sys.argv[1], sys.argv[2])

