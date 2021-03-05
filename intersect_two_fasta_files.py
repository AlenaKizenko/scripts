#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# choose all fasta ids, which are in first file, from second file
# arg1 = file with all repeats hmmer RT1, arg2 = second file with rest or expanded repeats
from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys


input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]

file2_set = {title for title, seq in SimpleFastaParser(open(input_file2))}

with open(output_file, 'w') as out_handle:
    with open(input_file1) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            if title in file2_set:
                out_handle.write(">%s\n%s\n" % (title, seq))