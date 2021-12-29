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
        te_class = line[2].split('_')[0]
        repid = re.search('(?<="Motif:).+?(?=")', line[8]).group()
        lst_rows.append(f'{line[0]}\t{line[1]}\tTE\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\tID={te_class}_{repid}\n')


with open(output_gff, 'w') as f_out:
    f_out.writelines(lst_rows) 