#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def split_seqid(cds, pep, output_folder):
    pep_lst = {}
    base_name = os.path.splitext(os.path.basename(os.path.normpath(pep)))[0]
    target_records = []
    for seq_record in SeqIO.parse(pep, "fasta"):
        type_rep = seq_record.id.split('_')[0]
        seq_record.id = seq_record.id.split('_')[1:]
        seq_record.id = '_'.join(seq_record.id)
        pep_lst[seq_record.id] = type_rep
    for seq_record2 in SeqIO.parse(cds, "fasta"):
        seq_record2.id = seq_record2.id.split('.p')[0]
        if seq_record2.id in pep_lst.keys():
            seq_record2.id = f'{pep_lst[seq_record2.id]}_{seq_record2.id}'
            rec = SeqRecord(seq_record2.seq, id=seq_record2.id, description='')
            target_records.append(rec)
    SeqIO.write(target_records, f'{output_folder}/{base_name}_cds.fasta', 'fasta')


result = split_seqid(sys.argv[1], sys.argv[2], sys.argv[3])