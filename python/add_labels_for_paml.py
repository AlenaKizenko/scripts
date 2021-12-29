#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from ete3 import Tree
import sys
import statistics
import regex as re


tree = Tree(sys.argv[1])
dist_lst = []

for node in tree.traverse():
    dist_lst.append(float(node.dist))
    me = statistics.mean(dist_lst)  # calculate distances mean
    lst_me = []  # create list for absolute deviations from mean
    for i in dist_lst:
        a = abs(i - me)  # from each distance subtract mean and take absolute value
        lst_me.append(a)  # append list with absolute deviation from mean
    mean_abs_dev = statistics.mean(lst_me) * 3  # return mean absolute deviation * 2

for node in tree.traverse():
    if node.is_leaf():
        try:
            re.search('hvir.*', node.name).group()
            node = node.up
            node.name = "#0"
        except:
            node.name = node.name
    else:
        node.name = ''

cnt = 1
for node in tree.traverse():
    if not node.is_leaf() and node.dist > mean_abs_dev and not node.name:
        node.name = f"#{cnt}"
        cnt += 1
    elif not node.is_leaf() and node.dist <= mean_abs_dev and not node.name:
        node.name = ""
    elif node.is_leaf():
        pass

for node in tree.traverse():
    if node.name:
        for i in node.children:
            if not i.name:
                i.name = node.name
                
for node in tree.traverse():                
    if not node.name:
        node.name = f'#{cnt+1}'

               
print(tree.write(format=8))