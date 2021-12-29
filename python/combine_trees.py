#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import os
import ete3


tree_omega = sys.argv[1]
tree_branch_length = sys.argv[2]

def create_dict_for_omega(tree_omega):
    tree_omega = ete3.Tree(tree_omega)
    #omega_dict = {}
    omega_lst = []
    for node in tree_omega.traverse():
        omega_lst.append(node.dist)
     #   omega_dict[node.name] = node.dist
    #return omega_dict
    return omega_lst

def add_omega(odict, tree_branch_length):
    cnt = 0
    tree_branch_length = ete3.Tree(tree_branch_length)
    for node in tree_branch_length.traverse():
        node.support = odict[cnt]
        cnt += 1
    return tree_branch_length




res_dict = create_dict_for_omega(sys.argv[1])
print(res_dict)
tree = add_omega(res_dict, sys.argv[2])
base_name = os.path.splitext(os.path.basename(os.path.normpath(sys.argv[2])))[0]
tree.write(format=0, outfile=f'{sys.argv[3]}/{base_name}_branch_omega.fasta')