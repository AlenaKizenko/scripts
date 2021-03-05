#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:04:30 2020

@author: alena
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import argcomplete
import os
#from BCBio import GFF
from itertools import islice
import re


def make_arguments_parser():
    parser = argparse.ArgumentParser(
        prog='Repeat organizer',
        description='Extracting particular repeat sequences from database and adding organism prefix')
        
    parser.add_argument('-ifa', '--input_fasta', 
                        help='Input FASTA file with repeats', 
                        required=True,
                        type=str)
    parser.add_argument('-igff', '--input_gff', 
                        help='Input GFF file with repeats', 
                        required=True,
                        type=str)                     
    parser.add_argument('-o','--output',
                        help='Output folder name', 
                        required=True,
                        type=str)
    argcomplete.autocomplete(parser) 
                                               
    return parser.parse_args()

if __name__ == "__main__":
    
    args = make_arguments_parser()
    
    try:
        os.mkdir(args.output)  # try to create folder
    except FileExistsError:
        print(f"{args.output} : directory already exists")  # if folder exists
    finally:
        os.chdir(args.output)  # change working directory to output folder
        path = os.getcwd()  # get working directory path
        print("Output directory path: {}".format(path))
    
    gff_ids = []
    
    with open(args.input_gff) as gff:
        for line in gff:
           # a = line.split('\t')[8]
            gff_ids.append(re.sub('ID=', '', line.split('\t')[8].split(';')[2]))
    
    renamed_records = []
    
    cnt = 0
    name = os.path.splitext(os.path.basename(os.path.normpath(args.input_fasta)))[0]
    for seq_record in SeqIO.parse(args.input_fasta, "fasta"):
        seq_record_id = f'{seq_record.id}_{gff_ids[cnt]}'
        cnt += 1
        rec = SeqRecord(seq_record.seq, id = seq_record_id, description = '')
        renamed_records.append(rec)
    SeqIO.write(renamed_records, f'{args.output}/{name}_fam.fasta', 'fasta')
