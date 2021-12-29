#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys
from Bio.SeqUtils import GC



input_fastq = sys.argv[1]
output_fastq = sys.argv[2]
#gc_cutoff = int(sys.argv[3])


with open(output_fastq, 'w') as out_handle:
    for title, seq, qual in FastqGeneralIterator(input_fastq):
        out_handle.write(f'{GC(seq)}\n')
        #if GC(seq) > gc_cutoff:
            #out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
