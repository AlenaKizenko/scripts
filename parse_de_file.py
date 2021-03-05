#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import re


gff = sys.argv[1]
genes_dict_condition = {}
genes_dict_de = {}

with open(gff) as file:
    file = file.readlines()
    for line in file[1:]:
        line = line.strip().split('\t')
        id_gene = re.search('(?<=ID=).+?(?=;)', line[8]).group()
        gene = f'{line[0]}.{id_gene}'
        genes_dict_condition[gene] = re.search('(?<=condition=).+?(?=;)', line[8]).group().split(',')
        genes_dict_de[gene] = re.search('(?<=de=).*', line[8]).group().split(',')

#print(genes_dict_de)

for key, value in genes_dict_condition.items():
    for count, i in enumerate(value):
        print(f'{key}\t{value[count]}\t{genes_dict_de[key][count]}')



            
    

