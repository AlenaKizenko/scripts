#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys


input_gff = sys.argv[1]
input_out = sys.argv[2]
output_gff = sys.argv[3]
lst_rows = []
names = {}


with open(input_out) as file:
    te_dict = {}
    file = file.readlines()
    for line in file[3::]:
        line = line.rstrip().split()
        te_dict[line[9]] = line[10]


with open(input_gff) as file:
    file = file.readlines()
    for line in file[3::]:
       # print(line)
        line = line.rstrip().split('\t')
        repid = re.search('(?<="Motif:).+?(?=")', line[8]).group()
        if not repid in names.keys():
            names[repid] = 1
        elif repid in names.keys():
            names[repid] += 1
        lst_rows.append(f'{line[0]}\t{line[1]}\t{te_dict[repid]}_{repid}_{names[repid]}\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\t{line[8]}\n')


with open(output_gff, 'w') as f_out:
    f_out.writelines(lst_rows) 