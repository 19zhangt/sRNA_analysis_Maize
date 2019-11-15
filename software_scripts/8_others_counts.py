#!/usr/bin/python
#-*-coding:utf-8 -*-

import os
import re
import sys
import copy
import subprocess
from PIL import Image

order = 0
pic_out = sys.argv[1]

show_value = {}
value_file = open(sys.argv[2], "r")
for seq_name in value_file.readlines():
    seq_name = seq_name.strip()
    seq_inform = seq_name.split("\t")
    # print(seq_inform)
    show_value[seq_inform[0]] = seq_inform[21]
value_file.close()

page_wid = 1600
page_hei = 1100

seq_str = {}
with open("data/6mireapout/forisoPremirnaSeq.fa") as pre_seq:
    for eli in pre_seq.readlines():
        if eli.startswith(">"):
            name = eli.strip()[1:]
            seq = ""
        else:
            seq += eli.strip()
            seq_str[name] = seq

for extract_name in list(show_value.keys())[order:]:
    over, seqs = subprocess.getstatusoutput('grep "' + extract_name + '\t" ' + sys.argv[3])
    if over == 0:
        seqs = seqs.split("\n")
        # print ("Read number:" + str(len(seqs)))
    seq_tpm = {}
    seq_loc = {}
    out_seq = []
    count_other = 0
    count_rev = 0
    for i in seqs:
        i_line = i.split("\t")
        if float(i_line[8]) > 0:
            seq_tpm[i_line[2]] = i_line[8]
            seq_loc[i_line[2]] = int(i_line[5])
            if i_line[1] == "others":
                count_other += 1
                if i_line[7] == "-":
                    count_rev += 1
    print(extract_name, show_value[extract_name], count_other, count_rev)
            # out_seq.append(i_line[1])
    # out_seq2 = sorted(seq_loc.items(), key=lambda item:item[1])
    # for i in out_seq2:
    #     out_seq.append(i[0])
    #     print(i)
