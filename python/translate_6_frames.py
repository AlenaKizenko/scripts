#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys


input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
aa_lst = []

for seq_record in SeqIO.parse(input_fasta, 'fasta'):
    dna_seqs = [seq_record.seq, seq_record.seq.reverse_complement()]
    aa_seqs = (s[i:].translate(to_stop=True) for i in range(3) for s in dna_seqs)
    for count, aa in enumerate(aa_seqs):
        aa_record = SeqRecord(aa, id = f'{seq_record.id}_{count}', description='')
        aa_lst.append(aa_record)

SeqIO.write(aa_lst, output_fasta, 'fasta')