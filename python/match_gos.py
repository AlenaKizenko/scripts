#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import csv


go_file = sys.argv[1]
res_file = sys.argv[2]
lst_gos = {}

with open(go_file) as file:
        file = file.readlines()
        for line in file:
            line = line.strip().split()
            if line[0] not in lst_gos.keys():
                lst_gos[line[0]] = []
            else:
                lst_gos[line[0]].append(line[1])


lst = []
for key, val in lst_gos.items():
    if len(val) != 0:
        val_pretty = ', '.join(val)
        lst.append(f'{key}\t{val_pretty}\n')

with open(res_file, 'w') as f_out:
    f_out.writelines(lst)