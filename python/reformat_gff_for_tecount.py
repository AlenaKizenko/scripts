#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys


input_gff = sys.argv[1]
output_gff = sys.argv[2]
lst_rows = []
names = {}



with open(input_gff) as file:
    file = file.readlines()
    for line in file:
       # print(line)
        line = line.rstrip().split('\t')
        gene_id = re.search('(?<="Motif:).+?(?=")', line[8]).group()
        class_family_id = line[2].split('_')[0]
        class_id = class_family_id.split('/')[0]
        try:
            family_id = class_family_id.split('/')[1]
        except:
            family_id = class_id
        num = line[2].split('_')[-1]
        transcript_id = f'{gene_id}_{num}'
        lst_rows.append(f'{line[0]}\t{line[1]}\texon\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\tgene_id {gene_id}; transcript_id {transcript_id}; family_id {family_id}; class_id {class_id};\n')


with open(output_gff, 'w') as f_out:
    f_out.writelines(lst_rows) 