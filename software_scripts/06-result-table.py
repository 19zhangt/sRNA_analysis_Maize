#! /usr/bin/env/ python
# _*_ coding:utf-8 _*_

import re
import os
import sys
import copy
import itertools
import subprocess


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
    return matchedPosList

def getMatchedPositions( startPos, endPos, matchedStructList ):
    RNALen = len(matchedStructList)
    #check pos
    # if( startPos < 0 or endPos < 0 or startPos > RNALen or endPos > RNALen ):
    #   print "Warnning: matrched pos not correctly determined:", startPos, endPos, RNALen
    #   return [-1,-1]
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
    if( matchedPosList[0] == -1 ):
        idx = 0
        for i in range(len(matchedPosList)):
            if( matchedPosList[i] != -1 ):
                idx = i
                break
        matchedPos = matchedPosList[idx]
        idxList = range(idx)
        for i in idxList[::-1]:
            matchedPos = matchedPos + 1
            matchedPosList[i] = matchedPos
    if( matchedPosList[-1] == -1 ):
        idxList = range(len(matchedPosList))
        idx = 0
        for i in idxList[::-1]:
            if( matchedPosList[i] != -1 ):
                idx = i
                break
        #end for i
        matchedPos =  matchedPosList[idx]
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
    if list(mirnaStruct).count("(") > len(mirnaStruct)/2:  # '(' in mirnaStruct:
        armDetailedType = "5p"
    if list(mirnaStruct).count(")") > len(mirnaStruct)/2:  # ')' in mirnaStruct:
        armDetailedType = "3p"
    return armDetailedType

def pre_mirna_to_table(name_ids, tab_pre_seq, tab_pre_str, tab_pre_loc, tab_mat_seq):
    """ The mature miRNA and miRNA* location information derived from pre-miRNA
    Returns:
        The table contains mature miRNA location information
    """
    rna_location = getMiRNAPosition(tab_mat_seq, tab_pre_seq)
    mirna_star = getMiRNAStar(tab_pre_seq, tab_pre_str, tab_mat_seq)
    rna_star_location = getMiRNAPosition(mirna_star, tab_pre_seq)
    tab_mat_seq_arm = checkArm(tab_pre_seq, tab_pre_str, tab_mat_seq)
    mirna_star_arm = checkArm(tab_pre_seq, tab_pre_str, mirna_star)
    pre_loc_info = re.split(":", tab_pre_loc)
    chromosome = pre_loc_info[0]
    strand = pre_loc_info[3]
    if strand == "+":
        mirna_location = '%s:%s:%s:%s' % (chromosome, str(int(pre_loc_info[1])+rna_location[0]),
                                          str(int(pre_loc_info[1])+rna_location[1]-1), strand)
        mirna_star_location = '%s:%s:%s:%s' % (chromosome, str(int(pre_loc_info[1])+rna_star_location[0]),
                                               str(int(pre_loc_info[1])+rna_star_location[1]-1), strand)
    else:
        mirna_location = '%s:%s:%s:%s' % (chromosome, str(int(pre_loc_info[2])-rna_location[1]+1),
                                          str(int(pre_loc_info[2])-rna_location[0]), strand)
        mirna_star_location = '%s:%s:%s:%s' % (chromosome, str(int(pre_loc_info[2])-rna_star_location[1]+1),
                                               str(int(pre_loc_info[2])-rna_star_location[0]), strand)
    output_result = ["miRName", tab_mat_seq, len(tab_mat_seq), tab_mat_seq_arm, mirna_location, "mat_abu", "conser", "miRSta", mirna_star, len(mirna_star), mirna_star_arm,
                     mirna_star_location, "mat_abu", "conser", name_ids]
    return output_result

""" start loading data
1. mireap-xxx.aln file contain the pre-miRNA name, location, sequence, length
"""
first_table = {}
with open(sys.argv[1]) as seq_len:
    for eli in seq_len.readlines():
        eli = eli.strip().split('\t')
        first_table[eli[0]] = eli
# aln_file = open(sys.argv[1], "r")
# for pre_name in aln_file.readlines():
#     line_split = re.split(" ", pre_name)
#     if len(line_split) >= 3:
#         nameInd = line_split[0]
#         nameInd2 = line_split[1]
#         if nameInd.startswith("xxx"):
#             first_table[nameInd] = [nameInd]
#             first_table[nameInd].append(line_split[1])
#         elif nameInd2 in first_table:
#             first_table[nameInd2].append(re.sub("T", "U", line_split[0]))
#             first_table[nameInd2].append(len(line_split[0]))
# aln_file.close()

structure_file = {}
raw = 1
structure_read = open(sys.argv[2], 'r')
for structure_line in structure_read.readlines():
    structure_line = structure_line.strip("\n")
    if structure_line.startswith(">"):
        structure_file[raw] = structure_line.strip(">")
    elif structure_line.startswith('.') or structure_line.startswith('('):
        structure_line = re.split(" ", structure_line)[0]
        structure_file[structure_file[raw]] = structure_line
        raw += 1
structure_read.close()


conInf = {}
second_table = {}
nameList = []
raw = 1
miran_star = open('data/9mirnas/ann_mirnas.fa', 'w')
pipeline_result = open(sys.argv[3], 'r')
for pre_name in pipeline_result.readlines():
    pre_name_split = pre_name.split("\t")
    if pre_name_split[0] == "Yes":
        name_id = pre_name_split[1]
        raw += 1
        nameList.append(name_id)
        mature_seq = pre_name_split[10].strip("\n")
        pre_seq = first_table[name_id][2]
        pre_location = first_table[name_id][1]
        if name_id not in second_table:
            new_string = pre_mirna_to_table(name_id, pre_seq, structure_file[name_id],pre_location, mature_seq)
            miran_star.write('>%s_%s\n%s\n>%s_%s\n%s\n' % (new_string[1], new_string[2],new_string[1],
                                                        new_string[8], new_string[9],new_string[8]))
            conInf[new_string[1]] = "No-conserved"
            conInf[new_string[8]] = "No-conserved"
            tmp_chr = first_table[name_id]
            tmp_chr.extend(["Abun", "Str_bias", "Abu_bias"])
            tmp_chr.extend(new_string)
            second_table[name_id] = first_table[name_id]
pipeline_result.close()
miran_star.close()
print(len(second_table), len(nameList))
"""
According to the table, find the known pre-miRNAs and mature miRNA and rename them
"""
known_pre = {}
# using the know name
currdir=os.path.dirname(sys.argv[7])
species= sys.argv[8]

known_mature = {}
mature_file = open('data/9mirnas/miRNAs.fa', 'r')
for mat in mature_file.readlines():
    mat = mat.strip()
    if mat.startswith(">"):
        pass
    else:
        known_mature[mat.replace("T", "U")] = "Know"
mature_file.close()


# The 14, 15, 16 name
nametpm = {}
resFil = open("data/9mirnas/getFinalRes_Cen.txt", "w")
nameover = {}
for allPre in second_table.keys():
    pre_name = second_table[allPre][0]
    miRseq = second_table[allPre][8]
    miRsta = second_table[allPre][15]
    nametpm[miRseq] = 0
    nametpm[miRsta] = 0
    if pre_name in known_pre:
        second_table[allPre][21] = known_pre[pre_name]
        if pre_name not in nameover:
            if miRseq in known_mature or miRsta in known_mature:
                second_table[allPre][7] = known_pre[pre_name]+"-"+second_table[allPre][10]
                second_table[allPre][14] = known_pre[pre_name]+"-"+second_table[allPre][17]
                nameover[pre_name] = 0
            else:
                second_table[allPre][7] = known_pre[pre_name]+".1-"+second_table[allPre][10]
                second_table[allPre][14] = known_pre[pre_name]+".1-"+second_table[allPre][17]
                nameover[pre_name] = 1
        else:
            if nameover[pre_name] == 0:
                second_table[allPre][7] = known_pre[pre_name]+".1-"+second_table[allPre][10]
                second_table[allPre][14] = known_pre[pre_name]+".1-"+second_table[allPre][17]
            if nameover[pre_name] == 1:
                second_table[allPre][7] = known_pre[pre_name]+".2-"+second_table[allPre][10]
                second_table[allPre][14] = known_pre[pre_name]+".2-"+second_table[allPre][17]
    myd = second_table[allPre]
    resFil.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\n'.
                 format(myd[0], myd[1], myd[2], myd[3], myd[4], myd[5], myd[6], myd[7], myd[8], myd[9], myd[10], myd[11],
                        myd[12], myd[13], myd[14], myd[15], myd[16], myd[17], myd[18], myd[19], myd[20], myd[21]))

resFil.close()

# mature miRNA coversation
abundance_file = open(sys.argv[4], "r")
for abu_line in abundance_file.readlines()[1::]:
    abu_line = abu_line.strip()
    abu_line_s = abu_line.split("\t")
    abu_line_s[2] = re.sub("T", "U", abu_line_s[2])
    nametpm[abu_line_s[2]] = round(float(abu_line_s[8]),2)
abundance_file.close()

# ssearch36 -m 8 -E 1  ./2-res/miRNA_input.fa  ./1-raw/miRBase22_miRNA.fa  >Result.txt
# mature miRNA conservation
# awk -F"\t" '{OFS="\t";split($2,c,"-");a=length($1);b=int($3*$4/100+0.5-$6);if(a-b<=2) print $1,c[1]}'  all_Result.txt | cut -f 1,2 | sort -k 1,1 -k 2,2 | uniq | awk -F"\t" '{a[$1]=a[$1]"""/"""$2}END{OFS="\t";for(i in a)print i,a[i]}' | awk -F"\t" 'NR==FNR{a[$1]=$2;next}{OFS="\t";if(a[$1]){print $1,a[$1]} else {print $1, "-"}}'  - all.txt |sed 's/\t\//\t/' >./2-res/Conserved1.txt
# ssearch36 -m 8 -E 1 ann_mirnas.fa mirbase_mature.fa  > Result.txt
# awk -F"\t" '{OFS="\t";split($1,c,"_");split($2,d,"-");a=c[2];b=int($3*$4/100+0.5-$6);if(a-b<=2) print c[1],d[1]}' Result.txt | sort -k 1,1 -k 2,2 | uniq | awk -F"\t" '{a[$1]=a[$1]"""/"""$2}END{OFS="\t";for(i in a)print i,a[i]}' | sed 's/\t\//\t/' >Conserved.txt
comm1 = 'ssearch36 -m 8 -E 1 data/9mirnas/ann_mirnas.fa data/9mirnas/miRNAs.fa >data/9mirnas/Result.txt'
a, b = subprocess.getstatusoutput(comm1)
comm2 = 'awk -F"\t" \'{OFS="\t";split($1,c,"_");split($2,d,"-");a=c[2];b=int($3*$4/100+0.5-$6);if(a-b<=2) print c[1],d[1]}\' data/9mirnas/Result.txt | sort -k 1,1 -k 2,2 | uniq | awk -F"\t" \'{a[$1]=a[$1]"""/"""$2}END{OFS="\t";for(i in a)print i,a[i]}\' | sed "s/\t\//\t/" >data/9mirnas/Conserved.txt'
a, b = subprocess.getstatusoutput(comm2)


con_res = open("data/9mirnas/Conserved.txt", "r")
for i in con_res.readlines():
    i = i.strip("\n")
    conseq = i.split("\t")[0]
    infrcon = i.split("\t")[1]
    infrcon = infrcon.split("/")
    infrcon = list(set(infrcon))
    infrcon = "/".join(infrcon)
    if infrcon == "-":
        conInf[conseq] = "No-conserved"
    else:
        conInf[conseq] = "Conserved"
con_res.close()


pre_abu_in = {}
pre_abu = open(sys.argv[5], "r")
for pre_line in pre_abu.readlines():
    pre_line = pre_line.strip()
    pre_info = pre_line.split("\t")
    pre_abu_in[pre_info[1]] = [pre_info[4], pre_info[2], pre_info[3]]
pre_abu.close()

if len(sys.argv)<10:
    finalName = {}
else:
    nameChange = open(sys.argv[9], "r")
    finalName = {}
    for i in nameChange.readlines():
        nameList = i.split("\t")
        nameList[1] = nameList[1].strip()
        finalName[nameList[0]] = nameList[1]
    nameChange.close()

print ('Write to '+ sys.argv[6] +'\n')
resFil = open(sys.argv[6], "w")

for allPre in second_table.keys():
    pre_name = second_table[allPre][0]
    miRseq = second_table[allPre][8]
    miRsta = second_table[allPre][15]
    second_table[allPre].append(nametpm[miRseq])
    second_table[allPre].append(nametpm[miRsta])
    second_table[allPre].append(conInf[miRseq])
    second_table[allPre].append(conInf[miRsta])
    if pre_name in finalName:
        second_table[allPre][21] = finalName[pre_name]
    else:
        finalName[pre_name] = pre_name
        second_table[allPre][21] = pre_name
    ## mature miRNAs name
    if finalName[pre_name].find("_N")!=-1:
        second_table[allPre][7] = finalName[pre_name]+"-"+second_table[allPre][10]
        second_table[allPre][14] = finalName[pre_name]+"-"+second_table[allPre][17]
    else:
        if miRseq in known_mature:
            second_table[allPre][7] = finalName[pre_name]+"-"+second_table[allPre][10]
        else:
            second_table[allPre][7] = finalName[pre_name]+".1-"+second_table[allPre][10]
        if miRsta in known_mature:
            second_table[allPre][14] = finalName[pre_name]+"-"+second_table[allPre][17]
        else:
            second_table[allPre][14] = finalName[pre_name]+".1-"+second_table[allPre][17]
    # else:
    #     if nameover[pre_name] == 0:
    #         second_table[allPre][7] = finalName[pre_name]+".1-"+second_table[allPre][10]
    #         second_table[allPre][14] = finalName[pre_name]+".1-"+second_table[allPre][17]
    #     if nameover[pre_name] == 1:
    #         second_table[allPre][7] = finalName[pre_name]+".2-"+second_table[allPre][10]
    #         second_table[allPre][14] = finalName[pre_name]+".2-"+second_table[allPre][17]
    myd = second_table[allPre]
    # information
    if myd[0] in pre_abu_in:
        # print myd[0]
        myd[12] = nametpm[myd[8]]
        myd[19] = nametpm[myd[15]]
        myd[13] = conInf[myd[8]]
        myd[20] = conInf[myd[15]]
        #
        myd[4] = round(float(pre_abu_in[myd[0]][0]),2)
        myd[5] = pre_abu_in[myd[0]][1]
        myd[6] = pre_abu_in[myd[0]][2]
        # print myd
        if myd[16]>18 and myd[16]<25:
            if(myd[9]>22 and myd[19]== 0):
                print(myd[0], myd[9])
                pass
            else:
                if(myd[12]>0 and myd[19]>0):
                    resFil.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\tHC\n'.
                                format(myd[0], myd[1], myd[2], myd[3], myd[4], myd[5], myd[6], myd[7], myd[8], myd[9], myd[10], myd[11],
                                        myd[12], myd[13], myd[14], myd[15], myd[16], myd[17], myd[18], myd[19], myd[20], myd[21]))
                else:
                    resFil.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\tLC\n'.
                                format(myd[0], myd[1], myd[2], myd[3], myd[4], myd[5], myd[6], myd[7], myd[8], myd[9], myd[10], myd[11],
                                        myd[12], myd[13], myd[14], myd[15], myd[16], myd[17], myd[18], myd[19], myd[20], myd[21]))

resFil.close()


# fas_structure = 'sort -k 9,9 ./0-FinalRes_Cen.txt| awk \'{print ">"""$1"""_"""$10"""\n"""$9}\'  >Tu_mature.fa'
# a, b = subprocess.getstatusoutput(fas_structure)
# ssearch36 -m 8 -E 1 Tu_mature.fa  mirbase_mature.fa >1result.txt
# awk -F"\t" '{OFS="\t";split($1,c,"_");split($2,d,"-");a=c[2];b=int($3*$4/100+0.5-$6);if(a-b<=2) print c[1],$2}' 1result.txt | sort -k 1,1 -k 2,2 | uniq | awk -F"\t" '{a[$1]=a[$1]"""/"""$2}END{OFS="\t";for(i in a)print i,a[i]}' | sed 's/\t\//\t/' >Conserved.txt
# fas_structure = 'sort -k 9,9 ./0-FinalRes_Cen.txt| awk \'{print ">"""$22"""_"""$10"""\n"""$9}\'  >Tu_mature_name2.fa'
# a, b = subprocess.getstatusoutput(fas_structure)
# ssearch36 -m 8 -E 1 Tu_mature_name2.fa Tu_mature_name2.fa > 2result.txt
# awk -F"\t" '{OFS="\t";split($1,c,"_");split($2,d,"_");a=c[2];b=int($3*$4/100+0.5-$6);if(a-b<=2&&c[1]!=d[1]) print c[1],d[1]}' 2result.txt | sort -k 1,1 -k 2,2 | uniq >Conserved3.txt
# grep ">" Tu_mature_name2.fa | sed -r 's/>//;s/_..//' | awk -F"\t" 'NR==FNR{a[$1]=$1;a[$2]=$2;next}{if(!a[$1]) print $1}' Conserved3.txt - >non_Conserved3.txt
