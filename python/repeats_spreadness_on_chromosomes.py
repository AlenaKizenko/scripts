#!/usr/bin/env python3
# -*- coding: utf-8 -*-



from Bio import SeqIO
import sys


genome = sys.argv[1]
gff_repeats = sys.argv[2]


def chr_length(genome):
    chr_len_dict = {}
    for seq_record in SeqIO.parse(genome, "fasta"):
        chr_len_dict[seq_record.id] = []
        chr_len_dict[seq_record.id].append(len(seq_record.seq))
    return chr_len_dict


def repeat_content(gff_repeats):
    gff_dict = {}
    with open(gff_repeats) as gff:
        gff = gff.readlines()
        for line in gff:
            line = line.strip().split()
            if not line[0] in gff_dict.keys():
                gff_dict[line[0]] = 1
            elif line[0] in gff_dict.keys():
                gff_dict[line[0]] += 1
    return gff_dict
    


chr_len = chr_length(genome)

repeats = repeat_content(gff_repeats)

for chrom in chr_len:
    if chrom in repeats.keys():
        chr_len[chrom].append(repeats[chrom])
    else:
        chr_len[chrom].append('0')

for key, value in chr_len.items():
    print(f'{key}\t{value[0]}\t{value[1]}')