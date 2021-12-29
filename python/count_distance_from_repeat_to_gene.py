#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import re


repeat_ann_file = sys.argv[1]

def filter_genes(file):
    with open(file) as file:
        file = file.readlines()
        for line in file:
            line = line.strip().split('\t')
            try:
                repid = re.search('(?<=Target=).+?(?= )', line[17]).group()
                print(line[17])
                if (int(line[3].strip()) - int(line[12].strip())) < 0: # repeat is upstream
                    distance = int(line[12].strip()) - int(line[4].strip())
                    if distance > 0:
                        print(f'{line[0]}\t{line[1]}\t{line[2]}\t'
                              f'{line[3]}\t{line[4]}\t{line[5]}\t'
                              f'{line[6]}\t{line[7]}\t{line[8]}\t'
                              f'{line[9]}\t{line[10]}\t{line[11]}\t'
                              f'{line[12]}\t{line[13]}\t{line[14]}\t'
                              f'{line[15]}\t{line[16]}\t{line[17]}\t'
                              f'{distance}\t{repid}')
                    elif distance < 0:
                        print(f'{line[0]}\t{line[1]}\t{line[2]}\t'
                              f'{line[3]}\t{line[4]}\t{line[5]}\t'
                              f'{line[6]}\t{line[7]}\t{line[8]}\t'
                              f'{line[9]}\t{line[10]}\t{line[11]}\t'
                              f'{line[12]}\t{line[13]}\t{line[14]}\t'
                              f'{line[15]}\t{line[16]}\t{line[17]}\t0\t{repid}')
                    else:
                        pass
                elif (int(line[3].strip()) - int(line[12].strip())) > 0: # repeat is downstream
                    distance = int(line[3].strip()) - int(line[13].strip())
                    if distance > 0:
                        print(f'{line[0]}\t{line[1]}\t{line[2]}\t'
                              f'{line[3]}\t{line[4]}\t{line[5]}\t'
                              f'{line[6]}\t{line[7]}\t{line[8]}\t'
                              f'{line[9]}\t{line[10]}\t{line[11]}\t'
                              f'{line[12]}\t{line[13]}\t{line[14]}\t'
                              f'{line[15]}\t{line[16]}\t{line[17]}\t'
                              f'-{distance}\t{repid}')
                    elif distance < 0:
                        print(f'{line[0]}\t{line[1]}\t{line[2]}\t'
                              f'{line[3]}\t{line[4]}\t{line[5]}\t'
                              f'{line[6]}\t{line[7]}\t{line[8]}\t'
                              f'{line[9]}\t{line[10]}\t{line[11]}\t'
                              f'{line[12]}\t{line[13]}\t{line[14]}\t'
                              f'{line[15]}\t{line[16]}\t{line[17]}\t0\t{repid}')
                    else:
                        pass
            except:
                pass

res = filter_genes(repeat_ann_file)