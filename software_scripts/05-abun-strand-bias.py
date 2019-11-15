
import re
import sys
import os.path
import subprocess
import pandas as pd

def isParenthese(ch):
    if( ch=='(' or ch==')' ):
        return True
    else:
        return False

def matchRNAStructure( RNAStructure ):
    RNALen = len(RNAStructure)
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

def get_seq_str(align_file):
    raw_file = open(align_file, "r")
    sequence_str = {}
    for align in raw_file.readlines():
        align_line = align.split(" ")
        if len(align_line) == 3:
            if re.match(r'xxx-m\d+$', align_line[1]):
                sequence = align_line[0]
                other, structure = subprocess.getstatusoutput("echo " + sequence + "| RNAfold --noPS")
                structure = re.split("\s", structure)[1]
                sequence_str[align_line[1]] = [sequence, structure]
    raw_file.close()
    return sequence_str


def strand_bias_iso(bias_name, value_list):
    iso_mi_rna = iso_result[bias_name]
    sense_val = 0
    sense_iso = 0
    nonsense = 0
    for sublist in value_list:
        sublist[6] = sublist[6].strip("\n")
        if sublist[4] == "+" and sublist[6] not in iso_mi_rna:
            sense_val += float(sublist[5])
    for iso in iso_mi_rna:
        if iso not in abundance:
            abundance[iso] = 0
        sense_iso += abundance[iso]
    for no_sublist in value_list:
        if no_sublist[4] == "-":
            nonsense += float(no_sublist[5])
    bias_val = (sense_val+sense_iso) / (sense_val+sense_iso+nonsense)
    return round(bias_val,3)


def abundance_bias_iso(bias_name, value_list):
    sense_list = []
    tpm = []
    for sublist in value_list:
        if sublist[4] == "+":
            sense_list.append(sublist)
            tpm.append(float(sublist[5]))
    # top 1 miRNA and his varieties
    iso_mi_rna = iso_result[bias_name]
    iso_abundance = 0
    for i in iso_mi_rna:
        iso_abundance += abundance[i]
    no_iso = 0
    for no_i in sense_list:
        seq_sub = no_i[6]
        seq_sub = seq_sub.strip("\n")
        if seq_sub not in iso_mi_rna:
            no_iso += float(no_i[5])
    bias_value = iso_abundance / (iso_abundance + no_iso)
    all_value = iso_abundance + no_iso
    return [iso_abundance, all_value, round(bias_value,3)]


def centroid_structure(rna_name, rna_seq):
    fas_file = open("tmp.fasta", "w")
    fas_file.write('>%s\n%s' % (rna_name, rna_seq))
    fas_file.close()
    fas_str = 'centroid_fold -e CONTRAfold -g 4 tmp.fasta > tmp.structure'
    other, res = subprocess.getstatusoutput(fas_str)
    if other == 0:
        str_file = open("tmp.structure", "r")
        for i in str_file.readlines():
            if i.startswith(".") or i.startswith("("):
                return i.split(" ")[0]
        str_file.close()

#re_str_seq, matNum
def stru_filter(test_seq, mature_mi_rna):
    mi_rna_pos = getMiRNAPosition(mature_mi_rna, test_seq[1])
    mi_star = getMiRNAStar(test_seq[1], test_seq[2], mature_mi_rna)
    mi_star_pos = getMiRNAPosition(mi_star, test_seq[1])
    pairStr = matchRNAStructure(test_seq[2])
    duplexS = mi_rna_pos[0]
    duplexE = mi_rna_pos[1]-2
    location=pairStr[duplexS:duplexE]
    location_bc=pairStr[duplexS:duplexE]
    print(location)
    while -1 in location_bc:
        location_bc.remove(-1)
    if all([location_bc[i] > location_bc[i+1] for i in range(len(location_bc)-1)]):
        if mi_rna_pos[0] > mi_star_pos[1]:
            loop_len = mi_rna_pos[0]-mi_star_pos[1]
            loop_ture = pairStr[(mi_star_pos[1]):(mi_star_pos[1]+3)]
        else:
            loop_len = mi_star_pos[0]-mi_rna_pos[1]
            loop_ture = pairStr[(mi_rna_pos[1]):(mi_rna_pos[1]+3)]
        print(test_seq[0],mi_rna_pos[0], mi_rna_pos[1], mi_star_pos[0], mi_star_pos[1], loop_len, loop_ture)
        if loop_len < 15 and loop_ture.count(-1) >0:
            return 'False' +'\t'+ str(loop_len) +'\t'+ str(loop_ture.count(-1)) +'\t'+ str('short_loop')
        elif(location.count(-1) >5):
            return 'False' +'\t'+ str('-') +'\t'+ str('-') +'\t'+ str('-')
        else:
            mismatch=0
            while location[0] == -1:
                mismatch+=1
                del location[0]
            bugleOne=0
            bugleTwo=0
            nomatch=0
            for num in range(len(location)-1):
                if(location.count(-1)>5):
                    mismatch = 99
                    break
                if(location[num]-location[num+1]==1):
                    nomatch=0
                elif(location[num]-location[num+1]==2):
                    print(num)
                    bugleOne+=1
                elif(location[num]-location[num+1]>2 and location[num] != -1 and location[num+1] != -1):
                    bugleTwo+=1
                    break
                elif(location[num+1] == -1):
                    nomatch+=1
                    if(location[num] != -1):
                        leftnum = location[num]
                    elif(num+1 == len(location)):
                        comstrand = leftnum-location[num+1]-1
                        mismatch = 99
                        break
                else:
                    comstrand = leftnum-location[num+1]-1
                    #print(location[num], nomatch, comstrand)
                    if(abs(nomatch-comstrand)==1):
                        bugleOne+=1
                    elif(abs(nomatch-comstrand)>1):
                        bugleTwo+=1
                    else:
                        mismatch+=min(nomatch, comstrand)
                        nomatch=0
            # print(mismatch,bugleOne,bugleTwo)
            if mismatch+bugleOne <= 5 and bugleOne <= 3 and bugleTwo == 0:
                return 'True' + '\t'+ str(mismatch+bugleOne) +'\t'+ str(bugleOne) +'\t'+ str(bugleTwo)
            elif(mismatch==99):
                return 'False' +'\t'+ str('-') +'\t'+ str('-') +'\t'+ str('-')
            else:
                return 'False' +'\t'+ str(mismatch+bugleOne) +'\t'+ str(bugleOne) +'\t'+ str(bugleTwo)
    else:
        print(location_bc)
        return 'False' +'\t'+ '-' +'\t'+ '-' +'\t'+ 'nonor'


## premiRNAs name and sequence
seq_str = {}
with open("data_miRNAs/6mireapout/forisoPremirnaSeq.fa") as pre_seq:
    for eli in pre_seq.readlines():
        if eli.startswith(">"):
            name = eli.strip()[1:]
            seq = ""
        else:
            seq += eli.strip()
            ot, structure = subprocess.getstatusoutput("echo " + seq + "| RNAfold --noPS")
            structure = re.split("\s", structure)[1]
            seq_str[name] = [seq, structure]
# seq_str = get_seq_str('data_miRNAs/6mireapout/mireap-xxx.aln')

## mature sequence
pre_to_mat = {}
# pre_mat = open('data_miRNAs/6mireapout/posMirnaSeq.txt', 'r')
# for mat_line in pre_mat.readlines():
#     mat_line = mat_line.strip()
#     mat_inform = mat_line.split("\t")
#     pre_to_mat[mat_inform[0]] = mat_inform[1]
#
# pre_mat.close()

print(seq_str[list(seq_str.keys())[0]])
bias_prename = []
out_data = {}
df_sran = pd.read_csv('data_miRNAs/7isomiRs/alignment_tpm_result.txt', sep='\t')
for pre_name in df_sran['pre_name'].unique():
    sub_frame = df_sran[df_sran['pre_name']==pre_name]
    if 'mature' in sub_frame['mir_type'].values:
        pre_to_mat[pre_name] = sub_frame[sub_frame['mir_type']=='mature']['sequence'].values[0]
        if sub_frame[sub_frame['strand']!='-'].shape[0] >1 and sum(sub_frame['RPM']>5)>1:
            ## strand bias
            print(sub_frame)
            mature_start = int(sub_frame[sub_frame['mir_type']=='mature']['start'])
            mature_end = int(sub_frame[sub_frame['mir_type']=='mature']['end'])
            bias_index = (abs(sub_frame['start'] - mature_start) <=3) & (abs(sub_frame['end'] - mature_end) <=3) & (sub_frame['strand'] == "-")
            mature_abu = int(sub_frame[sub_frame['mir_type']=='mature']['RPM'])
            strand_bias = round(mature_abu/(mature_abu+sum(sub_frame[bias_index]['RPM'])), 3)
            ## abundance bias
            sum_abun = sum(sub_frame['RPM'])
            true_abun = sum(sub_frame[sub_frame['strand']!='-']['RPM'])
            mat_index = (sub_frame['mir_type']!='others') & (sub_frame['strand']!='-')
            mat_iso = sum(sub_frame[mat_index]['RPM'])
            output = ['-',pre_name, str(strand_bias), str(round(mat_iso/true_abun,3)), str(mat_iso), str(true_abun)]
            if mat_iso/true_abun > 0.75 and strand_bias > 0.9:
                bias_prename.append(pre_name)
                output[0] = 'b'
                out_data[pre_name] = output
            elif sub_frame[sub_frame['mir_type']=='others'].shape[0] >1:
                min_frame = sub_frame[sub_frame['mir_type']=='others']
                second_seq = min_frame.loc[min_frame['RPM'].idxmax()]['sequence']
                second_star = getMiRNAStar(seq_str[pre_name][0], seq_str[pre_name][1], second_seq)
                if second_star in list(min_frame['sequence']):
                    start_mat = min_frame.loc[min_frame['RPM'].idxmax()]['start']
                    end_mat = min_frame.loc[min_frame['RPM'].idxmax()]['end']
                    con_one = (abs(min_frame['start'] - start_mat) <=3) & (abs(min_frame['end'] - end_mat) <=3)
                    start_star = min_frame[min_frame['sequence'] == second_star]['start']
                    end_star = min_frame[min_frame['sequence'] == second_star]['end']
                    con_two = (abs(min_frame['start'] - start_star) <=3) & (abs(min_frame['end'] - end_star) <=3)
                    add_value = sum(min_frame[con_one|con_two]['RPM'])
                    if (mat_iso+add_value)/true_abun > 0.75:
                        output[0] = 'ba'
                        output[3] = str(round((mat_iso+add_value)/true_abun,3))
                        print(pre_name)
                        out_data[pre_name] = output
                        bias_prename.append(pre_name)
                    else:
                        out_data[pre_name] = output
                else:
                    out_data[pre_name] = output
            else:
                out_data[pre_name] = output
        # else:
        #     out_data[pre_name] = ['-',pre_name, '-', '-', str(int(sub_frame[sub_frame['strand']!='-']['RPM'])),
        #                         str(int(sub_frame[sub_frame['strand']!='-']['RPM']))]


outfile = open('data_miRNAs/8bias/1-pipeline_res.txt', 'w')
for pre_name in bias_prename:
    re_str_vec = centroid_structure(pre_name, seq_str[pre_name][0])
    re_str_seq = [pre_name, seq_str[pre_name][0], re_str_vec]
    re_str_seq[1] = re_str_seq[1].replace("T", "U")
    mature_mi_rna_vec = pre_to_mat[pre_name]
    for matNum in [mature_mi_rna_vec]:
        matNum = matNum.replace("T", "U")
        # centroidFold filtering
        centroidfold_res = stru_filter(re_str_seq, matNum)
        centroidfold_tf = centroidfold_res.split("\t")[0]
        # RNAfold filtering
        rnafold_iuput=[re_str_seq[0], re_str_seq[1], seq_str[pre_name][1]]
        rnafold_res = stru_filter(rnafold_iuput, matNum)
        rnafold_tf = rnafold_res.split("\t")[0]
        # two tools result
        if rnafold_tf!= 'False' and centroidfold_tf != 'False' :
            out_data[pre_name][0] = out_data[pre_name][0] + '|R|C'
            outfile.write('Yes\t' + pre_name + '\t' + rnafold_res + '\t' + centroidfold_res + '\t' + matNum +'\n')
                # \t%s\t%s\t%s\t%s\n' % ([pre_name, rnafold_res, centroidfold_res, matNum]))
        elif rnafold_tf != 'False':
            out_data[pre_name][0] = out_data[pre_name][0] + '|R'
            outfile.write('No\t' + pre_name + '\t' + rnafold_res + '\t' + centroidfold_res + '\t' + matNum +'\n')
        elif centroidfold_tf != 'False':
            out_data[pre_name][0] = out_data[pre_name][0] + '|C'
            outfile.write('No\t' + pre_name + '\t' + rnafold_res + '\t' + centroidfold_res + '\t' + matNum +'\n')

outfile.close()


with open('data_miRNAs/8bias/2-bias_structure_info.txt', "w") as bias_pre:
    bias_pre.write('Type\tName\tStrand_bias\tAbundance_bias\tmiRNAs_isomiRs_RPM\tpre-miRNAs_RPM\n')
    for pre_name in out_data:
        bias_pre.write('\t'.join(out_data[pre_name]) + "\n")
