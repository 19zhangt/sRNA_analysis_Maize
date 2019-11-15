import re
import sys
import itertools

# mature and miRNA* location
location = {}
loc_file = open('data_miRNAs/6mireapout/posMirnaLocation.txt', 'r')
for loc_each in loc_file.readlines():
    loc_each = loc_each.strip()
    loc_line = loc_each.split('\t')
    location[loc_line[0]] = [int(loc_line[1]), int(loc_line[2]), int(loc_line[3]), int(loc_line[4])]

loc_file.close()

# results from isomiRs
def iso_type(ref_location,query_location):
    type_out = '-'
    if -3 <= query_location[0]-ref_location[0] < 0:
        type_out = 'add5'
        if query_location[1] > ref_location[1]:
            type_out = 'add5_add3'
        elif query_location[1] < ref_location[1]:
            type_out = 'add5_sub3'
    elif 0 < query_location[0]-ref_location[0] <= 3:
        type_out = 'sub5'
        if query_location[1] > ref_location[1]:
            type_out = 'sub5_add3'
        elif query_location[1] < ref_location[1]:
            type_out = 'sub5_sub3'
    elif 0 < query_location[1]-ref_location[1] <= 3:
        type_out = 'add3'
        if query_location[0] < ref_location[0]:
            type_out = 'add3_add5'
        elif query_location[0] > ref_location[0]:
            type_out = 'add3_sub5'
    elif -3 <= query_location[1]-ref_location[1] < 0:
        type_out = 'sub3'
        if query_location[0] < ref_location[0]:
            type_out = 'sub3_add5'
        elif query_location[0] > ref_location[0]:
            type_out = 'sub3_sub5'
    else:
        type_out = '-'
    return type_out

def iso_snp_type(ref_location,query_location,tmp_mutant):
    type_out = '-'
    seed_dist = int(tmp_mutant) + query_location[0] - ref_location[0]
    if query_location[1] - ref_location[1] == 1 and int(tmp_mutant) == query_location[1]-query_location[0] +1:
        type_out = 'nt_add3'
    if query_location[0] == ref_location[0] and query_location[1] == ref_location[1]:
        type_out = 'seed_snp'
    elif 9 <= seed_dist <= (ref_location[1]-ref_location[0]+1):
        type_out = 'tail_snp'
    return type_out

def iso_mat_star(mode, mat_location, star_location, query_location, mutant):
    type_overlap = ['-', '-']
    start_dis_1 = abs(query_location[0] - mat_location[0])
    end_dis_1 = abs(query_location[1] - mat_location[1])
    con_1 = (start_dis_1 <= 3 and end_dis_1 <= 3)
    start_dis_2 = abs(query_location[0] - star_location[0])
    end_dis_2 = abs(query_location[1] - star_location[1])
    con_2 = (start_dis_2 <= 3 and end_dis_2 <= 3)
    if mode == 'NonSNP':
        if start_dis_1 == 0 and end_dis_1 == 0:
            type_overlap = ['mature', '-']
        elif start_dis_1 <= 3 and end_dis_1 <= 3:
            type_overlap = ['iso',iso_type(mat_location, query_location)]
        elif start_dis_2 == 0 and end_dis_2 == 0:
            type_overlap = ['star', '-']
        elif start_dis_2 <= 3 and end_dis_2 <= 3:
            type_overlap = ['iso',iso_type(star_location, query_location)]
    else:
        snp_type_m = iso_snp_type(mat_location, query_location, mutant)
        snp_type_s = iso_snp_type(star_location, query_location, mutant)
        if start_dis_1 == 0 and end_dis_1 == 0 and snp_type_m == '-':
            type_overlap = ['mature', '-']
        elif start_dis_1 <= 3 and end_dis_1 <= 3 and snp_type_m != '-':
            type_overlap = ['iso', snp_type_m]
        elif start_dis_2 == 0 and end_dis_2 == 0 and snp_type_s == '-':
            type_overlap = ['star', '-']
        elif start_dis_2 <= 3 and end_dis_2 <= 3 and snp_type_s != '-':
            type_overlap = ['iso',snp_type_s]
    return type_overlap


mapping_result = {}
files_name = ['data_miRNAs/7isomiRs/Results/r0/RawResults/r0.txt', 'data_miRNAs/7isomiRs/Results/r1/RawResults/r1.txt']
for file_name in files_name:
    riso_file = open(file_name, "r")
    for iso_line in riso_file.readlines()[1::]:
        line_inform = iso_line.split("\t")
        line_new = line_inform[3]
        line_new = re.sub("\[\'", "", line_new)
        line_new = re.sub("\'\]", "", line_new)
        line_inform[3] = line_new.split("\', \'")
        pre_length = len(line_inform[3])
        iso_seq = line_inform[0]
        for i in range(pre_length):
            pre_each = line_inform[3][i].split(":")
            pre_name = pre_each[0]
            if pre_name in location:
                start_site = int(pre_each[1])
                end_site = start_site + int(len(line_inform[0])-1)
                if file_name == files_name[0]:
                    type_out_name = iso_mat_star('NonSNP', location[pre_name][0:2], location[pre_name][2:4], [start_site,end_site], '-')
                    if type_out_name[0] != '-':
                        if pre_name not in mapping_result:
                            mapping_result[pre_name] = [":".join([pre_name,type_out_name[0],iso_seq,type_out_name[1],'-',str(start_site),str(end_site),'+'])]
                            # all_alignment[pre_name] = [iso_seq + "\t" + str(int(pre_each[1]) -1)]
                        else:
                            mapping_result[pre_name].append(":".join([pre_name,type_out_name[0],iso_seq,type_out_name[1],'-',str(start_site),str(end_site),'+']))
                            # all_alignment[pre_name].append(iso_seq + "\t" + str(int(pre_each[1]) -1))
                else:
                    mutant_type = pre_each[2].split(",")[0]
                    mutant = pre_each[2].split(",")[1]
                    type_out_name = iso_mat_star('SNP', location[pre_name][0:2], location[pre_name][2:4], [start_site,end_site], mutant)
                    if type_out_name[0] != '-':
                        if pre_name not in mapping_result:
                            mapping_result[pre_name] = [":".join([pre_name,type_out_name[0],iso_seq,type_out_name[1],mutant_type,str(start_site),str(end_site),'+'])]
                        else:
                            mapping_result[pre_name].append(":".join([pre_name,type_out_name[0],iso_seq,type_out_name[1],mutant_type,str(start_site),str(end_site),'+']))
    riso_file.close()

# abundance_seq = list(itertools.chain.from_iterable(mapping_result.values()))
# abundance_seq = list(set(abundance_seq))
# iso_file = open('0-unique-isomiRNAs.txt', "w")
# for i in abundance_seq:
#     iso_file.write(i + "\n")
# iso_file.close()
# isomiR_file = open('pre-miRNAs2isomiRs.txt', "w")
# for i in mapping_result:
#     for each in list(set(mapping_result[i])):
#         isomiR_file.write('%s\t%s\n' % (i, each))
# isomiR_file.close()

mapping_file = open('data_miRNAs/7isomiRs/0-alignment.txt', "r")
for aline in mapping_file:
    line_detail = aline.strip().split('\t')
    pre_name = line_detail[1]
    match_list = [pre_name,'others',line_detail[6],'-','-',line_detail[2],line_detail[3],line_detail[4]]
    if pre_name not in mapping_result:
        mapping_result[pre_name] = [":".join(match_list)]
    elif line_detail[6] not in [x.split(':')[2] for x in mapping_result[pre_name]]:
        mapping_result[pre_name].append(":".join(match_list))
    else:
        pass

mapping_file.close()

# using plot the show
alignment_file = open('data_miRNAs/7isomiRs/0-pre-miRNAs_mapping_read.txt', "w")
for i in mapping_result.values():
    for each_loc in i:
        alignment_file.write('%s\n' % each_loc)

alignment_file.close()
