#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re

input_repeat_file = sys.argv[1]
input_deg_file = sys.argv[2]
deg_lst = []


with open(input_deg_file) as file:
    file = file.readlines()
    for line in file[1:]:
        line = line.strip().split('\t')
        id_gene = re.search('(?<=ID=).+?(?=;)', line[8]).group()
        gene = f'{line[0]}.{id_gene}'
        deg_lst.append(gene)


with open(input_repeat_file) as file:
        file = file.readlines()
        for line in file:
            line = line.strip().split('\t')
            #print(line[17])
            try:
                repid = re.search('(?<=Target=).+?(?= )', line[17]).group()
                if repid in deg_lst:
                    print(f'{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t'\
                             f'{line[7]}\t{line[8]}\t{line[9]}\t{line[10]}\t{line[11]}\t{line[12]}\t{line[13]}\t'\
                             f'{line[14]}\t{line[15]}\t{line[16]}\t{line[17]}\t{repid}')
            except:
                pass