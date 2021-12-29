#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import ete3
import sys


def label_nodes(tree):
    tree = ete3.Tree(tree)
    cnt = 1
    for leaf in tree.iter_leaves():  # iterate on leaves
        if not leaf.name:
            leaf.name = cnt
            cnt += 1
    print(tree)
    tree.write(format=1, outfile= '/Users/alenakizenko/Documents/PhD/hydra_vulgaris_chrom/RT1_analysis/tree.nw')

res = label_nodes(sys.argv[1])
