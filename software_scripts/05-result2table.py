#! /usr/bin/env/ python
# _*_ coding:utf-8 _*_

import re
import os
import sys
import subprocess

""" start loading data
1. mireap-xxx.aln file contain the pre-miRNA name, location, sequence, length
"""

seq_str = {}
with open("data/6mireapout/forisoPremirnaSeq.fa") as pre_seq:
    for eli in pre_seq.readlines():
        if eli.startswith(">"):
            name = eli.strip()[1:]
            seq = ""
        else:
            seq += eli.strip()
            seq_str[name] = seq

"""Generating the pre-miRNA structure file
Using the aln file and pipeline_result to extract pre-miRNA structure
"""
tmpfile='data/8bias/pipe_pre-miRNA.fas'
tmp_name = []
tmp_file = {}
pre_fas_file = open(tmpfile, 'w')
pipeline_result = open('data/8bias/1-pipeline_result.txt', 'r')
for pre_name in pipeline_result.readlines():
    pre_name_split = pre_name.split('\t')
    if pre_name_split[0] == "Yes":
        name_id = pre_name_split[1]
        tmp_name.append(name_id)
        pre_seq = seq_str[name_id]
        if name_id not in tmp_file:
            tmp_file[name_id] = 1
            pre_fas_file.write('>%s\n%s\n' % (name_id, pre_seq))
        else:
            pass

pipeline_result.close()
pre_fas_file.close()


## Sequence structures
fas_structure = 'RNAfold --noPS --noLP -i '+tmpfile+' > %s' % ('data/8bias/3-pipe_pre-miRNA.str')
a, b = subprocess.getstatusoutput(fas_structure)
fas_structure = 'sed -i -r "s/ .+//" %s' % ('data/8bias/3-pipe_pre-miRNA.str')
a, b = subprocess.getstatusoutput(fas_structure)


first_table = {}
aln_file = open('data/6mireapout/mireap-xxx.aln', "r")
for pre_name in aln_file.readlines():
    line_split = re.split(" ", pre_name)
    if len(line_split) >= 3:
        nameInd = line_split[0]
        nameInd2 = line_split[1]
        if nameInd.startswith("xxx"):
            first_table[nameInd] = [nameInd]
            first_table[nameInd].append(line_split[1])
        elif nameInd2 in first_table:
            first_table[nameInd2].append(re.sub("T", "U", line_split[0]))
            first_table[nameInd2].append(len(line_split[0]))
aln_file.close()

with open("genome/miRBasemiRNAs/02_miRNA_seq_len.txt") as seq_len:
    for eli in seq_len.readlines():
        eli = eli.strip().split('\t')
        eli[2] = re.sub("T", "U", eli[2])
        first_table[eli[0]] = eli

with open('data/9mirnas/mir_seq_len.txt', 'w') as write_seq:
    for key in first_table:
        if key in tmp_name:
            out_tmp = [str(x) for x in first_table[key]]
            write_seq.write("\t".join(out_tmp) + "\n")
