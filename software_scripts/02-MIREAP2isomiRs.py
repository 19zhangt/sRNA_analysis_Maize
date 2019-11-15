#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract most abundant tags as candidate for each precursor
Author: Ting Zhang
"""

import re
import subprocess
import os
import sys

def isParenthese(ch):
    if( ch=='(' or ch==')' ):
        return True
    else:
        return False

def matchRNAStructure( RNAStructure ):
    RNALen = len(RNAStructure)
    #first define a list with the same length of RNAStructure, at each position
    #the value of matched position for current position
    matchedPosList = [-1]*RNALen
    stack = []
    stackPos = []
    for i in range(RNALen):
        a = RNAStructure[i]
        if( isParenthese(a) == False ):
            continue
        if( len(stack) == 0 ):
            #empty stack, record item
            stack.append(a)
            stackPos.append(i)
        else:
            #pop last item in stack and stackPos
            stack_lastItem = stack.pop()
            stackPos_lastItem = stackPos.pop()
            if( stack_lastItem == '(' and a == ')' ):
                #meet a match, record matched position
                matchedPosList[i] = stackPos_lastItem
                matchedPosList[stackPos_lastItem] = i
                continue
            else:
                #have to record
                stack.append( stack_lastItem )
                stackPos.append( stackPos_lastItem )
                stack.append(a)
                stackPos.append(i)
            #end else
        #end else
    #end for i
    return matchedPosList

def getMatchedPositions( startPos, endPos, matchedStructList ):
    RNALen = len(matchedStructList)
    #check pos
    # if( startPos < 0 or endPos < 0 or startPos > RNALen or endPos > RNALen ):
    #     print "Warnning: matrched pos not correctly determined:", startPos, endPos, RNALen
    #     return [-1,-1]
    #end if
    if( startPos < 0 ):
        startPos = 0
    if( endPos < 0 ):
        endPos = endPos
    if( startPos > RNALen ):
        startPos = RNALen - 1
    if( endPos > RNALen ):
        endPos = RNALen - 1
    if( startPos == 0 and endPos == 0 ):
        return [0,0]
    matchedPosList = []
    for i in range(startPos, endPos):
        curPos = matchedStructList[i]
        matchedPosList.append( curPos )
    #end for i
    if( matchedPosList[0] == -1 ):
        idx = 0
        for i in range(len(matchedPosList)):
            if( matchedPosList[i] != -1 ):
                idx = i
                break
        #end for i
        matchedPos = matchedPosList[idx]
        idxList = range(idx)
        for i in idxList[::-1]:
            matchedPos = matchedPos + 1
            matchedPosList[i] = matchedPos
        #end for i
    #end if
    if( matchedPosList[-1] == -1 ):
        idxList = range(len(matchedPosList))
        idx = 0
        for i in idxList[::-1]:
            if( matchedPosList[i] != -1 ):
                idx = i
                break
        #end for i
        matchedPos =    matchedPosList[idx]
        idx = idx + 1
        for i in range(idx, len(matchedPosList)):
            matchedPos = matchedPos - 1
            matchedPosList[i] = matchedPos
    #end if
    if( matchedPosList[0] < 0 ):
        matchedPosList[0] = 0
    if( matchedPosList[0] > RNALen ):
        matchedPosList[0] = RNALen - 1
    if( matchedPosList[-1] < 0 ):
        matchedPosList[-1] = 0
    if( matchedPosList[-1] > RNALen ):
        matchedPosList[-1] = RNALen - 1
    return [matchedPosList[0], matchedPosList[-1]]

def getMiRNALen (miRNASeq):
    return len(miRNASeq)

def getMiRNAPosition (miRNASeq, RNASequence):
    startPos = RNASequence.find(miRNASeq)
    endPos = startPos + getMiRNALen(miRNASeq)
    return [startPos, endPos]

def getMiRNAStar(RNASequence, RNAStructure, miRNASeq, matchedStructList = [] ):
    if( len(matchedStructList) == 0 ):
        matchedStructList = matchRNAStructure( RNAStructure )
    [miRNAPosStart, miRNAPosEnd] = getMiRNAPosition( miRNASeq, RNASequence )
    [matchedPosStart, matchedPosEnd] = getMatchedPositions( miRNAPosStart, miRNAPosEnd-2, matchedStructList )
    matchedPosStart = matchedPosStart + 2
    matchedPosEnd = matchedPosEnd
    RNALen = len(RNASequence)
    if( matchedPosStart < 0 ):
        matchedPosStart = 0
    if( matchedPosEnd < 0 ):
        matchedPosEnd = 0
    if( matchedPosStart > RNALen ):
        matchedPosStart = RNALen - 1
    if( matchedPosEnd > RNALen ):
        matchedPosEnd = RNALen - 1
    starSeq = ""
    if( matchedPosStart < matchedPosEnd ):
        starSeq = RNASequence[matchedPosStart: (matchedPosEnd+1)]
        starSeq = starSeq[::-1]
    else:
        starSeq = RNASequence[matchedPosEnd: (matchedPosStart+1)]
    return starSeq

def getMiRNAStructure (miRNASeq, RNASequence, RNAStructure ):
    [startPos, endPos] = getMiRNAPosition( miRNASeq, RNASequence )
    mirnaStruct = RNAStructure[ startPos:endPos ]
    return mirnaStruct

def getMirnaIncludedInLoop( RNASequence, RNAStructure, miRNASeq ):
    flag = False
    mirnaStruct = getMiRNAStructure( miRNASeq, RNASequence, RNAStructure )
    if( ('(' in mirnaStruct) and (')' in mirnaStruct) ):
        flag = True
    return flag

def checkArm( RNASequence, RNAStructure, miRNASeq ):
    armDetailedType = ""
    if getMirnaIncludedInLoop(RNASequence, RNAStructure, miRNASeq):
        armDetailedType = "loop"
    #check arms
    mirnaStruct = getMiRNAStructure(miRNASeq, RNASequence, RNAStructure)
    armDetailedType = "unmatchedRegion"
    if list(mirnaStruct).count("(") > len(mirnaStruct)/2:    # '(' in mirnaStruct:
        armDetailedType = "5p"
    if list(mirnaStruct).count(")") > len(mirnaStruct)/2:    # ')' in mirnaStruct:
        armDetailedType = "3p"
    return armDetailedType

def uni_reap_gff(gff_input_file, gff_output_file):
    remain_pre=[]
    raw_len_range = range(0,0)
    raw_len=0
    line_strand=0
    raw_file = open(gff_input_file, "r")
    for gff_line in raw_file.readlines():
        gff_line_inform = gff_line.split("\t")
        if gff_line_inform[2] == "precursor":
            gff_pre_name = re.search('xxx-m\d+', gff_line_inform[8]).group()
            dist_len_range = range(int(gff_line_inform[3]), int(gff_line_inform[4]) + 1)
            dist_len=int(gff_line_inform[4])-int(gff_line_inform[3])+1
            overlap_len = list(set(dist_len_range).intersection(set(raw_len_range)))
            if len(overlap_len) > 1 and gff_line_inform[6]==line_strand:
                if dist_len>=raw_len:
                    remain_pre[len(remain_pre)-1]=gff_pre_name
                    line_start = int(gff_line_inform[3])
                    line_end = int(gff_line_inform[4])
                    raw_len_range = range(line_start, line_end)
                    line_strand=gff_line_inform[6]
                    raw_len=line_end-line_start+1
            else:
                remain_pre.append(gff_pre_name)
                line_start = int(gff_line_inform[3])
                line_end = int(gff_line_inform[4])
                raw_len_range = range(line_start, line_end)
                line_strand=gff_line_inform[6]
                raw_len=line_end-line_start+1
    raw_file.close()

    res_file = open(gff_output_file, "w")
    raw_file = open(gff_input_file, "r")
    for gff_line in raw_file.readlines():
        gff_line_inform = gff_line.split("\t")
        gff_pre_name = re.search('xxx-m\d+', gff_line_inform[8]).group()
        if gff_pre_name in remain_pre:
            res_file.write(gff_line)

    raw_file.close()
    res_file.close()
    return remain_pre


def get_seq_str(align_file):
    raw_file = open(align_file, "r")
    sequence_str = {}
    for align in raw_file.readlines():
        align_line = align.split(" ")
        if len(align_line) == 3:
            if re.match(r'xxx-m\d+$', align_line[1]):
                sequence = align_line[0]
                other, structure = subprocess.getstatusoutput("echo " + sequence + "| RNAfold")
                structure = re.split("\s", structure)[1]
                sequence_str[align_line[1]] = [sequence, structure]
    raw_file.close()
    return sequence_str


def get_mirna_star(star_pre_seq, star_mature_seq):
    a, b = subprocess.getstatusoutput('echo ' + star_pre_seq + '|RNAfold --noPS')
    if a == 0:
        pre_structure = re.split("\s", b)[1]
        mirna_star = getMiRNAStar(star_pre_seq, pre_structure, star_mature_seq)
        return mirna_star


"""
1. What we have files is the aln and gff files
remove redundance in gff file
"""
work_dir = sys.argv[1]
aln_file = work_dir + '/mireap-xxx.aln'
gff_file = work_dir + '/mireap-xxx.gff'
gff_unique = work_dir + '/mireap-xxx_uniq.gff'
pre_name_uni = uni_reap_gff(gff_file, gff_unique)

#1.1 The preparation of isomiRs files

iso_premirna = open(work_dir + '/forisoPremirnaSeq.fa', "w")
iso_mature_miRNA = open(work_dir + '/forisoMirnaSeq.fa', "w")

location = {}
mature_seq = {}
pre_to_mat = {}
raw = 0
max_name = 0
max_name_re = 0
max_seq = 0
max_star = 0
pre_name_back = 0
pre_seq_back = 0
max_value = 0
pre_name = 0

rawFile = open(aln_file, "r")
for i in rawFile.readlines():
    splList = i.strip().split(" ")
    if len(splList) == 3:
        if re.match(r'xxx-m\d+$', splList[1]):
            mature_seq[max_name] = max_seq
            if splList[1] in pre_name_uni:
                raw = 0
                pre_name = splList[1]
                pre_seq = splList[0]
                max_value = 0
                print(max_seq, len(str(max_star)), pre_seq_back, max_name, max_name_re, pre_name_back)
                if len(str(max_seq)) <= 26 and len(str(max_star)) <= 26 and pre_seq_back != 0 and max_name in ['3p', '5p'] and max_name_re in ['3p', '5p'] and pre_name_back not in pre_to_mat:
                    iso_mature_miRNA.write('>%s\n%s\n>%s\n%s\n' % (pre_name_back + "-" + max_name,max_seq, pre_name_back + "-" + max_name_re, max_star))
                    iso_premirna.write('>%s\n%s\n' % (pre_name_back, pre_seq_back))
                    print(max_seq, pre_seq_back)
                    mi_rna_l1 = getMiRNAPosition(max_seq, pre_seq_back)
                    mi_rna_l2 = getMiRNAPosition(max_star, pre_seq_back)
                    location[pre_name_back] = [mi_rna_l1[0], mi_rna_l1[1], mi_rna_l2[0], mi_rna_l2[1]]
                    pre_to_mat[pre_name_back] = max_seq
            else:
                    raw = 1
        elif re.match(r'mireap_', splList[1]) and raw == 0:
            print (splList, max_value, max_seq)
            if int(splList[2]) >= max_value:
                a, b = subprocess.getstatusoutput('echo ' + pre_seq + '| RNAfold --noPS')
                if a == 0:
                    pre_structure = re.split("\s", b)[1]
                    max_value = int(splList[2])
                    max_seq = splList[0].strip("-")
                    # pre_structure = centroid_structure(pre_name, pre_seq)
                    max_star = getMiRNAStar(pre_seq, pre_structure, max_seq)
                    # max_star = get_mirna_star(pre_seq, max_seq)
                    max_name = checkArm(pre_seq, pre_structure, max_seq)
                    max_name_re = checkArm(pre_seq, pre_structure, max_star)
                    pre_name_back = pre_name
                    pre_seq_back = pre_seq
                    #print(max_seq,pre_seq_back)

#print(max_seq,pre_seq_back)
mi_rna_l1 = getMiRNAPosition(max_seq, pre_seq_back)
mi_rna_l2 = getMiRNAPosition(max_star, pre_seq_back)
location[pre_name_back] = [mi_rna_l1[0], mi_rna_l1[1], mi_rna_l2[0], mi_rna_l2[1]]
pre_to_mat[pre_name_back] = max_seq
iso_mature_miRNA.write('>%s\n%s\n>%s\n%s\n' % (pre_name_back + "-" + max_name, max_seq, pre_name_back + "-" + max_name_re, max_star))
iso_premirna.write('>%s\n%s\n' % (pre_name_back, pre_seq_back))

rawFile.close()
iso_mature_miRNA.close()
iso_premirna.close()

## Saving the location
mature_loc = open(work_dir + '/posMirnaLocation.txt', "w")
for each_loc in location:
    res_tmp = location[each_loc]
    mature_loc.write('%s\t%s\t%s\t%s\t%s\n' % (each_loc, res_tmp[0]+1,res_tmp[1],res_tmp[2]+1,res_tmp[3]))
mature_loc.close()

## Saving the most abundance miRNA
mature_seq = open(work_dir + '/posMirnaSeq.txt', "w")
for mir in pre_to_mat:
    mature_seq.write('%s\t%s\n' % (mir, pre_to_mat[mir]))
mature_seq.close()
