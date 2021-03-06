#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


input_gff = sys.argv[1]
output_gff = sys.argv[2]
lst_names = sys.argv[3]
lst_rows = []
names = []
test = []

with open(lst_names) as clusters:
    clusters = clusters.readlines()
    for line in clusters:
        line = line.rstrip('\n').lstrip('Hvulgaris_')
        names.append(line)


with open(input_gff) as file:
    file = file.readlines()
    for line in file:
        line = line.strip().split('\t')
        if line[11].strip() in names:
            lst_rows.append(f'{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t'\
                             f'{line[7]}\t{line[8]}\t{line[9]}\t{line[10]}\t{line[11]}\t{line[12]}\t{line[13]}\t'\
                             f'{line[14]}\t{line[15]}\t{line[16]}\t{line[17]}\t{line[18]}\n')
       
       
with open(output_gff, 'w') as f_out:
    f_out.writelines(lst_rows) 