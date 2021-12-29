#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

# to get one largest contig 

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
contig_len = 0
contig_seq = str()
contig_name = str()
contig_des = str()


for seq_record in SeqIO.parse(input_fasta, "fasta"):
    if len(seq_record.seq) > contig_len:
        contig_len = len(seq_record.seq)
        contig_seq = seq_record.seq
        contig_name = seq_record.id
        contig_des = seq_record.description

rec = SeqRecord(contig_seq, id=contig_name, description=contig_des)
SeqIO.write(rec, output_fasta, 'fasta')

# to get 3 largest contigs

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
contigs = {}


for seq_record in SeqIO.parse(input_fasta, "fasta"):
    contigs[seq_record.id] = []
    contigs[seq_record.id].append(len(seq_record.seq))
    contigs[seq_record.id].append(seq_record.seq)
    contigs[seq_record.id].append(seq_record.description)
    
lst_len = []
for key, value in contigs.items():
    lst_len.append(value[0])

lst_len = sorted(lst_len, reverse=True)

lst_hits = []

for key, value in contigs.items():
    if value[0] == lst_len[0] or value[0] == lst_len[3] or value[0] == lst_len[2]:
        rec = SeqRecord(value[1], id=key, description=value[2])
        lst_hits.append(rec)

SeqIO.write(lst_hits, output_fasta, 'fasta')