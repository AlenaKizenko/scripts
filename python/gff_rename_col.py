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
    for line in file[2::]:
       # print(line)
        line = line.rstrip().split('\t')
        repid = re.search('(?<=;ID=).+?(?=;)', line[8]).group()
        if not repid in names.keys():
            names[repid] = 1
        elif repid in names.keys():
            names[repid] += 1
        lst_rows.append(f'{line[0]}\t{line[1]}\t{line[2]}_{repid}_{names[repid]}\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\t{line[8]}\n')


with open(output_gff, 'w') as f_out:
    f_out.writelines(lst_rows) 