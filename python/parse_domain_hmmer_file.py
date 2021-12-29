#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os


def cds_from_hmmmer(hmmer_file, cds_file, dom_file, output_folder):
    base_name = os.path.splitext(os.path.basename(os.path.normpath(cds_file)))[0]
    hmmer_res = SearchIO.read(hmmer_file, 'hmmsearch3-domtab')
    hmmer_dict = {}
    target_records = []
    for i in range(len(hmmer_res)-1):
        id_res = hmmer_res[i].id
        res = hmmer_res[i][0][0]
        hmmer_dict[id_res] = []
        hmmer_dict[id_res].append(res.hit_start)
        hmmer_dict[id_res].append(res.hit_end)
    
    dom_ids = [seq_record.id.split('.p')[0] for seq_record in SeqIO.parse(dom_file, 'fasta')]
    for seq_record in SeqIO.parse(cds_file, "fasta"):
        if seq_record.id in hmmer_dict.keys() and seq_record.id.split('.p')[0] in dom_ids:
            rec = SeqRecord(seq_record.seq[hmmer_dict[seq_record.id][0]*3:hmmer_dict[seq_record.id][1]*3], id=seq_record.id.split('.p')[0], description='')
            target_records.append(rec)        
    SeqIO.write(target_records, f'{output_folder}/{base_name}_domain.fasta', 'fasta')
    
lres = cds_from_hmmmer(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])