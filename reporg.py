#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import argcomplete
import os


def make_arguments_parser():
    parser = argparse.ArgumentParser(
        prog='Repeat organizer',
        description='Extracting particular repeat sequences from database and adding organism prefix')
        
    parser.add_argument('-i', '--input', 
                        help='Input FASTA file with repeats', 
                        required=True,
                        type=str)
    parser.add_argument('-pre', '--prefix',
                        help='Prefix (genome) name for adding at the beginning of each sequence',
                        required=True,
                        type=str)                     
    parser.add_argument('-o','--output',
                        help='Output folder name', 
                        required=True,
                        type=str)
    parser.add_argument('-rep', '--repeat',
                        required = True, type=str,
                        help='Target repeat name for filtering')
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
    
    target_records = []
    
    if '/' in args.repeat:
        rep = args.repeat.split('/')
        rep = '_'.join(rep)
    else:
        rep = args.repeat
    
    for seq_record in SeqIO.parse(args.input, "fasta"):
        if args.repeat in seq_record.id:
            seq_record_id = f'{args.prefix}_{seq_record.id}'
            rec = SeqRecord(seq_record.seq, id = seq_record_id, description = '')
            target_records.append(rec)
    SeqIO.write(target_records, f'{args.output}/{args.prefix}_{rep}.fasta', 'fasta')
