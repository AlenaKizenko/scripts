#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import re


gene_ann_file = sys.argv[1]
#dist_num = sys.argv[2]
wnt_genes = ['Sc4wPfr_1061.g18842.t1', 'Sc4wPfr_134.g20117.t1', 'Sc4wPfr_14.g1752.t1', 'Sc4wPfr_172.1.g1591.t1', 'Sc4wPfr_179.g7820.t1',\
             'Sc4wPfr_287.2.g28262.t1', 'Sc4wPfr_319.g27364.t1', 'Sc4wPfr_399.g28064.t1', 'Sc4wPfr_440.g6790.t1', 'Sc4wPfr_624.g29781.t1',\
             'Sc4wPfr_834.1.g27891.t1', 'Sc4wPfr_850.1.g5726.t1', 'Sc4wPfr_975.g7262.t1']


#wnt_genes = ['Sc4wPfr_1061.g18842.t1', 'Sc4wPfr_134.g20117.t1', 'Sc4wPfr_179.g7820.t1',\
 #            'Sc4wPfr_440.g6790.t1', 'Sc4wPfr_834.1.g27891.t1',' Sc4wPfr_850.1.g5726.t1'] 
def filter_genes(file):
    with open(file) as file:
        file = file.readlines()
        for line in file:
            line = line.strip().split('\t')
            repid = re.search('(?<=Target=).+?(?= )', line[8]).group()
            if repid in wnt_genes:
                print(f'{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\t{line[8]}\t')


res = filter_genes(gene_ann_file)          