#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import ete3


def classify_by_species(tree):
    for leaf in tree.iter_leaves():
        if 'Hvulgaris' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'Hydra_vulgaris'
                leaf = leaf.up
        elif 'Hviridissima' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'Hydra_viridissima'
                leaf = leaf.up
        elif 'Holigactis' in str(leaf):
            while leaf.up:
                if not str(leaf.name):
                    leaf.name = 'Hydra_oligactis'
                leaf = leaf.up
    return tree


tree = ete3.Tree('/home/alena/Documents/PhD/H_vir_oli_vulg_tree')
#print(tree)
name_species = classify_by_species(tree)
name_species.write(format=1, outfile="/home/alena/Documents/PhD/H_vir_oli_vulg_tree_named_nodes")
#print(name_species)
