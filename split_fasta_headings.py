#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# split fasta ids in hmm results by .p (Transdecoder numeration)
# arg1 = input filename, arg2 = output path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os


def split_seqid(filename, output_folder):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(filename)))[0]
    target_records = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq_record.id = seq_record.id.split(' ')[0]
        rec = SeqRecord(seq_record.seq, id=seq_record.id, description='')
        target_records.append(rec)
    SeqIO.write(target_records, f'{output_folder}/{base_name}_renamed.fasta', 'fasta')

result = split_seqid(sys.argv[1], sys.argv[2])