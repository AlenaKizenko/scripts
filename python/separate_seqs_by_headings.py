#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# split fasta ids in hmm results by .p (Transdecoder numeration)
# arg1 = input filename, arg2 = output path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import re


def split_seqid(filename, output_folder):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(filename)))[0]
    rest_records = []
    expansion_records = []
    for seq_record in SeqIO.parse(filename, "fasta"):
        if re.match('rest', seq_record.id):
            seq_record.id = f'{seq_record.id[14:]}_rest'
            rec = SeqRecord(seq_record.seq, id=seq_record.id, description='')
            rest_records.append(rec)
        elif re.match('expansion', seq_record.id):
            seq_record.id = f'{seq_record.id[19:]}_exp'
            rec = SeqRecord(seq_record.seq, id=seq_record.id, description='')
            expansion_records.append(rec)
    SeqIO.write(expansion_records, f'{output_folder}/{base_name}_expansion.fasta', 'fasta')
    SeqIO.write(rest_records, f'{output_folder}/{base_name}_rest.fasta', 'fasta')

result = split_seqid(sys.argv[1], sys.argv[2])