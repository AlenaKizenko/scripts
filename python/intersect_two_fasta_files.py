#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# choose all fasta ids, which are in first file, from second file
# arg1 = file with domains hmmer RT1, arg2 = file with repeats
from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys


input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]

file2_set = {title.split('.')[0] for title, seq in SimpleFastaParser(open(input_file2))}

with open(output_file, 'w') as out_handle:
    with open(input_file1) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            if title.split('.')[0] in file2_set:
                out_handle.write(">%s\n%s\n" % (title.split('.')[0], seq))