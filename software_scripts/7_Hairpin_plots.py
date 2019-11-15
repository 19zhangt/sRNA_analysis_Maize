#!/usr/bin/python
#-*-coding:utf-8 -*-

import os
import re
import sys
import copy
import subprocess
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.pdfbase.pdfmetrics import stringWidth

def format_values_ps(value):
    val_tab = [" %s 0 0.8 0 cfmark" % v for i, v in enumerate(value)]
    return "".join(val_tab)

def text_location(text, font, size, color, width, height):
    t = pdf.beginText()
    t.setFont(font, size)
    t.setCharSpace(0)
    t.setFillColor(color)
    t.setTextOrigin(width, height)
    t.textLine(text)
    return t

def text_line(text, width, color):
    t = pdf.beginText()
    t.setFont('Courier', 8)
    t.setCharSpace(0)
    t.setFillColor(color)
    t.setTextOrigin(width, raw+20)
    student_name = text
    text_width = stringWidth(student_name, 'Courier', 8)
    t.textLine(student_name)
    return t, text_width

name_change, show_value, iso_info = {}, {}, {}
with open("data_miRNAs/10plot/miRNAs_and_locations.txt") as locat:
    for eli in locat.readlines():
        eli = eli.strip().split("\t")
        name_change[eli[1]] = eli[0]

with open("data_miRNAs/10plot/finalTable.txt") as info_table:
    for eli in info_table.readlines():
        eli = eli.strip().split("\t")
        if eli[1] in name_change:
            real_name = name_change[eli[1]]
            name_change[eli[0]] = real_name
            show_value[real_name] = eli

with open("data_miRNAs/10plot/alignment_tpm_result.txt") as iso_mat:
    for eli in iso_mat.readlines():
        eli = eli.strip().split("\t")
        if eli[0] in name_change and eli[1]!="others":
            mirna_name = name_change[eli[0]]
            if mirna_name not in iso_info:
                iso_info[mirna_name] = [eli]
            else:
                iso_info[mirna_name].append(eli)

pic_out = "data_miRNAs/10plot/"
page_wid = 1600
page_hei = 1100
order = 1

for each_name in list(show_value.keys()):
    sequence = show_value[each_name][2]
    # Sequence to layout
    if len(sequence) < 200:
        seq_site = 900 - (stringWidth(sequence, 'Courier', 8)/2)
        raw = page_hei-100
        iso_number = 95
    else:
        seq_site = 770 - (stringWidth(sequence, 'Courier', 8)/2)
        raw = page_hei/2+50
        iso_number = 55

    # RNAfold struture
    file_name = pic_out+each_name
    other, structure = subprocess.getstatusoutput("echo " + sequence + "| RNAfold --noPS")
    structure = re.split("\s", structure)
    mfe_val = re.sub("\(|\)","",structure[2])
    inputFile = file_name+'.str'
    with open(inputFile, "w") as fas_file:
        fas_file.write(">"+each_name+"\n")
        fas_file.write(sequence+"\n")
        fas_file.write(structure[1] + " " + structure[2])

    # sequence start site
    tab_info = show_value[each_name]
    pre_loc = tab_info[1].split(':')
    mature_loc = tab_info[11].split(':')
    star_loc = tab_info[18].split(':')
    if pre_loc[3] == '+':
        mat_loc = [int(mature_loc[1])-int(pre_loc[1])+1, int(mature_loc[2])-int(pre_loc[1])+1,
        int(star_loc[1])-int(pre_loc[1])+1, int(star_loc[2])-int(pre_loc[1])+1]
    else:
        mat_loc = [int(pre_loc[2])-int(mature_loc[2])+1, int(pre_loc[2])-int(mature_loc[1])+1,
        int(pre_loc[2])-int(star_loc[2])+1, int(pre_loc[2])-int(star_loc[1])+1]
    # site in mature and miRNA*
    if tab_info[10] == "5p":
        sign1_name = "Mature"
        sign2_name = "Star"
        sign1_col = "red"
        sign2_col = "green"
        sign1_rgd = ['0.8', '0', '0']
        sign2_rgd = ['0', '0.8', '0']
        sign1 = int((mat_loc[0] + mat_loc[1])/2)
        sign2 = int((mat_loc[2] + mat_loc[3])/2)
        sign1_loc = seq_site + stringWidth(sequence[:sign1], 'Courier', 8) - stringWidth(sign1_name, 'Courier', 14)/2
        sign2_loc = seq_site + stringWidth(sequence[:sign2], 'Courier', 8) - stringWidth(sign2_name, 'Courier', 14)/2
    else:
        sign1_name = "Star"
        sign2_name = "Mature"
        sign1_col = "green"
        sign2_col = "red"
        sign1_rgd = ['0', '0.8', '0']
        sign2_rgd = ['0.8', '0', '0']
        sign2 = int((mat_loc[2] + mat_loc[3])/2)
        sign1 = int((mat_loc[0] + mat_loc[1])/2)
        sign1_loc = seq_site + stringWidth(sequence[:sign2], 'Courier', 8) - stringWidth(sign2_name, 'Courier', 14)/2
        sign2_loc = seq_site + stringWidth(sequence[:sign1], 'Courier', 8) - stringWidth(sign1_name, 'Courier', 14)/2
    # color
    mat_sort_loc = sorted(mat_loc)
    color_tab = []
    for i in range(len(sequence)):
        if mat_sort_loc[0] <= i+1 <= mat_sort_loc[1]:
            color_tab.append(str(i+1))
            color_tab.extend(sign1_rgd)
            color_tab.append('cfmark')
        elif mat_sort_loc[2] <= i+1 <= mat_sort_loc[3]:
            color_tab.append(str(i+1))
            color_tab.extend(sign2_rgd)
            color_tab.append('cfmark')
        elif mat_sort_loc[1] < i+1 < mat_sort_loc[2]:
            color_tab.extend([str(i+1), '1', '0.6', '0', 'cfmark'])
        else:
            color_tab.extend([str(i+1), '0.8', '0.8', '0.8', 'cfmark'])
    extraMacro = "/cfmark {setrgbcolor newpath 1 sub coor exch get aload pop fsize 2 div 0 360 arc fill} bind def"

    os.system("RNAplot --pre \" %s %s\" < %s" % (extraMacro, " ".join(color_tab), inputFile))
    os.system('convert -trim -quality 500  -density 500 ' + each_name + '_ss.ps ' + file_name + '.png && rm ' + each_name + '_ss.ps ')

    img = Image.open( file_name + '.png')
    (x, y) = img.size
    scale = float(max(x, y))/min(x, y)

    pdf = canvas.Canvas(pic_out+ str(order) + ".pdf", pagesize=(page_wid, page_hei))
    model = sequence
    model = sequence+"    "+str(round(float(tab_info[4]), 3))
    mat_sort_loc.insert(0, 0)
    mat_sort_loc.append(len(model))
    color_list = ["black", sign1_col, "#CD950C", sign2_col, "black"]
    raw_width = seq_site

    for i in range(len(mat_sort_loc)-1):
        if (i % 2) == 0:
            start = mat_sort_loc[i]
            end = mat_sort_loc[i+1]-1
        else:
            start = mat_sort_loc[i]-1
            end = mat_sort_loc[i+1]
        subseq = model[start:end]
        t, textWidth = text_line(subseq, raw_width, color_list[i])
        raw_width += textWidth
        pdf.drawText(t)

    t = pdf.beginText()
    t.setFont('Courier', 14)
    t.setCharSpace(0)
    t.setFillColor(sign1_col)
    t.setTextOrigin(sign1_loc, raw+30)
    t.textLine(sign1_name)

    t.setFillColor(sign2_col)
    t.setTextOrigin(sign2_loc, raw+30)
    t.textLine(sign2_name)

    t.setFont('Courier', 8)
    t.setFillColor("black")
    t.setTextOrigin(seq_site, raw+10)
    t.textLine(structure[1])

    t.setFillColor("red")
    t.setTextOrigin(seq_site, raw)
    mat_infor = [x for x in iso_info[each_name] if x[1]=='mature'][0]
    text_left = int(mat_infor[5])-1
    text_right = len(sequence)-int(mat_infor[6])
    only_seq = "."*text_left + re.sub("T","U",mat_infor[2]) + "."*text_right
    t.textLine(only_seq+"    "+mat_infor[8])
    star_infor = [x for x in iso_info[each_name] if x[1]=='star']
    if len(star_infor) == 1:
        raw -= 10
        t.setFillColor("green")
        t.setTextOrigin(seq_site, raw)
        text_left = int(star_infor[0][5])-1
        text_right = len(sequence)-int(star_infor[0][6])
        only_seq = "."*text_left + re.sub("T","U",star_infor[0][2]) + "."*text_right
        t.textLine(only_seq+"    "+star_infor[0][8])

    iso_seq = [x for x in iso_info[each_name] if x[1]=='iso']
    if len(iso_seq) > 0:
        t.setFillColor("black")
        if len(iso_seq) > iso_number:
            tmp_val = sorted([float(x[8]) for x in iso_seq], reverse= True)[iso_number]
            iso_seq = [x for x in iso_seq if float(x[8])>= tmp_val]
        if tab_info[10] == "5p":
            iso_seq = sorted(iso_seq, key=lambda x:int(x[5]))
        else:
            iso_seq = sorted(iso_seq, key=lambda x:int(x[5]), reverse= True)
        for eli in iso_seq:
            # if float(eli[8]) >= 10:
            raw -= 10
            t.setTextOrigin(seq_site, raw)
            text_left = int(eli[5])-1
            text_right = len(sequence)-int(eli[6])
            only_seq = "."*text_left + re.sub("T","U",eli[2]) + "."*text_right
            t.textLine(only_seq+"    "+eli[8])

    pdf.drawText(t)

    if len(sequence) < 200:
        pdf.drawImage(file_name + '.png', 150, 200, width=(400/scale), height=400)
        word1_site = 100
        word1_height = page_hei-150
    else:
        pdf.drawImage(file_name + '.png', 1200-(300/scale), page_hei-350, width=(300/scale), height=300)
        word1_site = 300
        word1_height = page_hei-100

    loc_info = tab_info[1].split(":")
    title = text_location(each_name, 'Courier', 15, "black", word1_site, word1_height)
    chr = text_location("Chromosome:    chr" + loc_info[0], 'Courier', 10, "black", word1_site, word1_height-40)
    sta = text_location("Start:         " + loc_info[1], 'Courier', 10, "black", word1_site, word1_height-50)
    end = text_location("End:           " + loc_info[2], 'Courier', 10, "black", word1_site, word1_height-60)
    strand = text_location("Strand:        " + loc_info[3], 'Courier', 10, "black", word1_site, word1_height-70)
    sbias = text_location("Strand bias:   " + str(round(float(tab_info[5]), 3)),
                            'Courier', 10, "black", word1_site, word1_height-80)
    abias = text_location("Abundance bias:" + str(round(float(tab_info[6]), 3)),
                            'Courier', 10, "black", word1_site, word1_height-90)
    mfe = text_location("MFE:           " + mfe_val, 'Courier', 10, "black", word1_site, word1_height-100)
    pdf.drawText(title)
    pdf.drawText(chr)
    pdf.drawText(sta)
    pdf.drawText(end)
    pdf.drawText(strand)
    pdf.drawText(sbias)
    pdf.drawText(abias)
    pdf.drawText(mfe)
    pdf.showPage()
    pdf.save()
    order += 1
