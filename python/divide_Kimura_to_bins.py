#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
import sys
import regex as re
import os


bin1, bin2, bin3, bin4 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
repeat_fasta = sys.argv[5]
output_folder = sys.argv[6]
output_folder = os.path.abspath(output_folder)


def read_bin_files(bin_file):
    bin_lst = []
    with open(bin_file) as file:
        file = file.readlines()
        for line in file:
            line = line.strip()
            bin_lst.append(line)
    return bin_lst

bin1_lst = read_bin_files(bin1)
bin2_lst = read_bin_files(bin2)
bin3_lst = read_bin_files(bin3)
bin4_lst = read_bin_files(bin4)

result_bin1, result_bin2, result_bin3, result_bin4 = str(), str(), str(), str()
result_records = []

for seq_record in SeqIO.parse(repeat_fasta, "fasta"):
    if re.search('(?<=hvul_LINE/CR1_).+?(?=_\w+\()', seq_record.id).group() in bin1_lst:
        result_bin1+=seq_record.seq
    elif re.search('(?<=hvul_LINE/CR1_).+?(?=_\w+\()', seq_record.id).group() in bin2_lst:
        result_bin2+=seq_record.seq
    elif re.search('(?<=hvul_LINE/CR1_).+?(?=_\w+\()', seq_record.id).group() in bin3_lst:
        result_bin3+=seq_record.seq
    elif re.search('(?<=hvul_LINE/CR1_).+?(?=_\w+\()', seq_record.id).group() in bin4_lst:
        result_bin4+=seq_record.seq

rec_bin1 = SeqRecord(result_bin1, id = 'LINE_CR1_bin1', description = '')
result_records.append(rec_bin1)
rec_bin2 = SeqRecord(result_bin2, id = 'LINE_CR1_bin2', description = '')
result_records.append(rec_bin2)
rec_bin3 = SeqRecord(result_bin3, id = 'LINE_CR1_bin3', description = '')
result_records.append(rec_bin3)
rec_bin4 = SeqRecord(result_bin4, id = 'LINE_CR1_bin4', description = '')
result_records.append(rec_bin4)

SeqIO.write(result_records, f'{output_folder}/LINE_CR1_binned.fasta', 'fasta')


