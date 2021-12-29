#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:15:35 2020

@author: alena
"""

import sys

def parse_blast_output(path:str):
    d = {}
    with open(path,"r") as f:
        for line in f:
            [query,subject,_,_,_,_,_,_,_,_,_,_,dist] = line.strip().split("\t")
            # concat query and subject (defined order)
            if query < subject:
                key = query + ";" + subject
            else:
                key = subject + ";" + query
            # add/update dict
            if key not in d:
                d[key] = dist

    # print
    print("Source\tTarget\tWeight")
    for k,v in d.items():
        ids = k.split(";")
        print(f"{ids[0]}\t{ids[1]}\t{v}")          


if __name__ == "__main__":
    parse_blast_output(sys.argv[1])

