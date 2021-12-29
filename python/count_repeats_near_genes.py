#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import csv

repeat_ann_file = sys.argv[1]
csv_file = sys.argv[2]
gene_dict = {}


with open(repeat_ann_file) as file:
    file = file.readlines()
    for line in file:
        line = line.strip().split('\t')
        if line[19] in gene_dict.keys():
            gene_dict[line[19]] += 1
        else:
            gene_dict[line[19]] = 1

with open(csv_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        for key, value in gene_dict.items():
            writer.writerow([key, value])           