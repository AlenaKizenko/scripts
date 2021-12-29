#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import os
import regex as re

in_dir = sys.argv[1]
output_file = sys.argv[2]
out_file = []


for file in os.listdir(in_dir):
    d = os.path.join(in_dir, file)
    if os.path.isdir(d):
        folder = os.path.basename(os.path.normpath(d)).split('~')[0]
        model_type, branch = folder.split('.')[0], folder.split('.')[1]
        with open(f'{d}/out') as out:
            out = out.readlines()
            for line in out:
                if line.startswith('lnL'):
                    np = re.search('(?<=np:).+?(?=\):)', line.strip()).group()
                    lnL = re.search('(?<=\):[\s]+).+?(?=[\s]+\+)', line.strip()).group()
                elif line.startswith('w (dN/dS)'):
                    w = line.strip().split()[5]
                    out_file.append(f'{model_type}\t{branch}\t{lnL}\t{np}\t{w}\n')

with open(output_file, 'w') as f_out:
    f_out.writelines(out_file) 