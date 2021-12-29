#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import os
import ete3


def split_seqid(new, old, tree):
    tree = ete3.Tree(tree)
    with open(new, 'r') as new:
        new_lst = [line.strip() for line in new.readlines()]
    with open(old, 'r') as old:
        old_lst = [line.strip() for line in old.readlines()]
    for leaf in tree.iter_leaves():
        if leaf.name in new_lst or leaf.name[3::] in new_lst:
            leaf.name = f'new_{leaf.name}'
        elif leaf.name in old_lst or leaf.name[3::] in old_lst:
            leaf.name = f'old_{leaf.name}'
    return tree


result = split_seqid(sys.argv[1], sys.argv[2], sys.argv[3])
base_name = os.path.splitext(os.path.basename(os.path.normpath(sys.argv[3])))[0]
result.write(format=1, outfile=f'{sys.argv[4]}/{base_name}_prefix.fasta')