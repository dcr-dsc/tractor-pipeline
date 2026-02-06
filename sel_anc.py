#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import argparse
import ast
import os

parser = argparse.ArgumentParser()
parser.add_argument('-ancs', type=str, required=True, help='Insert Ancestries')
# example: -ancs "'AFR','EUR'"
args = parser.parse_args()


try:
    sel_anc = ast.literal_eval(f"[{args.ancs}]")
except (SyntaxError, ValueError):
    print("Invalid ancestry format.")
    exit(1)

# Load metadata file with ancestry information
sample_map = pd.read_csv("b38_sample_map.txt", sep='\t', header=None)

sample_map.columns=['id', 'anc']

# Selecting Ancestries 
sel_code='.'.join([s.lower() for s in sel_anc])

sel_map=sample_map[sample_map['anc'].isin(sel_anc)]

sel_map.to_csv('b38_%s_map.txt' %(sel_code),sep='\t',header=False,index=False)
