#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def split_seqid(dom, pep, output_folder, type_seq):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(dom)))[0]
    domain_lst = []
    target_records = []
    for seq_record in SeqIO.parse(dom, "fasta"):
        domain_lst.append(seq_record.id.split('.p')[0])
    for seq_record2 in SeqIO.parse(pep, "fasta"):
        seq_record2.id = seq_record2.id.split('.p')[0]
        if seq_record2.id in domain_lst:
            rec = SeqRecord(seq_record2.seq, id=seq_record2.id, description='')
            target_records.append(rec)
    SeqIO.write(target_records, f'{output_folder}/{base_name}_{type_seq}.fasta', 'fasta')


result = split_seqid(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])