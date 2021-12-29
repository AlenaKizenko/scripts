#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

genome = sys.argv[1]
output_name = sys.argv[2]
target_records = []

for seq_record in SeqIO.parse(genome, "fasta"):
    seq_record.id = seq_record.id.split(' ')[0]
    rec = SeqRecord(seq_record.seq, id = seq_record.id, description = '')
    target_records.append(rec)

SeqIO.write(target_records, output_name, 'fasta')