
# mkdir data_miRNAs/7isomiRs data_miRNAs/8bias data_miRNAs/9mirnas
#
# awk -v refdir=data_miRNAs/6mireapout/forisoPremirnaSeq.fa -v matdir=data_miRNAs/6mireapout/forisoMirnaSeq.fa \
# '{if($0~/^lib:/){print "lib:", "data_miRNAs/5largematrix/reads_18_26.fa", "Sample1", "fa"} else if($0~/^main_ref:/) \
# {print "main_ref:", refdir, "yes"} else if($0~/^known_miRNAs:/) \
# { print "known_miRNAs:", matdir} else{ print $0}}' software_scripts/3_Config_bak.txt >data_miRNAs/7isomiRs/Config.txt
#
# python2 software_scripts/3_isomiRs_ID.py data_miRNAs/7isomiRs
#
#
# bowtie-build data_miRNAs/6mireapout/forisoPremirnaSeq.fa data_miRNAs/6mireapout/forisoPremirnaSeq
#
# bowtie -p 10 -v 0 -f -t -a data_miRNAs/6mireapout/forisoPremirnaSeq data_miRNAs/5largematrix/reads_18_26.fa >data_miRNAs/7isomiRs/Remapping.sam
#
# awk -v RS=">" 'NR>1{OFS="\t";print $1,$3}' data_miRNAs/5largematrix/reads_18_26.fa >data_miRNAs/5largematrix/reads_18_26.txt
#
# awk  'NR==FNR{a[$1]=$2;next}{OFS="\t"; print $1,$4,$5+1,$5+length($6),$3,$2,a[$1]}' \
# data_miRNAs/5largematrix/reads_18_26.txt data_miRNAs/7isomiRs/Remapping.sam | sort -k 2,2  >data_miRNAs/7isomiRs/0-alignment.txt
#
# python software_scripts/04-isomiRs-abundance.py
#
# sed -i "s/:/\t/g" data_miRNAs/7isomiRs/0-pre-miRNAs_mapping_read.txt
# awk -v RS=">" 'NR>1{OFS="\t";print $1,$2,$3}' data_miRNAs/5largematrix/reads_18_26.fa >data_miRNAs/7isomiRs/Uniq_Seq_TPM_sum.txt
#
# awk 'NR==FNR{a[$3]=$3;next}{OFS="\t";if(a[$3]){print $3,$2}}' \
# data_miRNAs/7isomiRs/0-pre-miRNAs_mapping_read.txt data_miRNAs/7isomiRs/Uniq_Seq_TPM_sum.txt >data_miRNAs/7isomiRs/tmp_alignment_read.txt

# awk 'BEGIN{OFS="\t";print "pre_name","mir_type","sequence","iso_type","mutant","start","end","strand","RPM"}NR==FNR{a[$1]=$2;next}{if(a[$3]){print $0"\t"a[$3]}}' \
# data_miRNAs/7isomiRs/tmp_alignment_read.txt data_miRNAs/7isomiRs/0-pre-miRNAs_mapping_read.txt > data_miRNAs/7isomiRs/alignment_tpm_result.txt

# python software_scripts/05-abun-strand-bias.py

# (head -n 1 data_miRNAs/8bias/2-bias_structure_info.txt && tail -n +2 data_miRNAs/8bias/2-bias_structure_info.txt | sort -k 1,1r ) > data_miRNAs/8bias/2-bias_structure_result.txt

# sort -k 1,1r -k 3,3r -k 7,7r data_miRNAs/8bias/1-pipeline_res.txt | awk 'BEGIN{OFS="\t";print "True","Name","RNAfold_true","RNAfold_mis","RNAfold_bugle",
# "RNAfold_abbugle","centeriodfold_true","centeriodfold_mis","centeriodfold_bugle","centeriodfold_abbugle","matureSeq"}{OFS="\t";print $0}' -  >data_miRNAs/8bias/1-pipeline_result.txt

# python software_scripts/05-result2table.py

# grep "miRNA_primary_transcript" genome/miRBasemiRNAs/zma.gff3 | sed 's/chr//' >data_miRNAs/9mirnas/Pre-miRNA.gff3
# awk '{print ">"$1"\n"$2}' genome/miRBasemiRNAs/mature_rmzma.fa >data_miRNAs/9mirnas/miRNAs.fa

# python software_scripts/06-result-table.py data_miRNAs/9mirnas/mir_seq_len.txt data_miRNAs/8bias/3-pipe_pre-miRNA.str \
# data_miRNAs/8bias/1-pipeline_result.txt data_miRNAs/7isomiRs/alignment_tpm_result.txt data_miRNAs/8bias/2-bias_structure_result.txt  \
# data_miRNAs/9mirnas/tmpTable1.txt data_miRNAs/9mirnas/ Zea_mays;

# awk -F"\t" '{OFS="\t";split($2, a,":");print a[1],a[2]-1,a[3],$1,".",a[4]}' data_miRNAs/9mirnas/tmpTable1.txt \
#  | bedtools intersect -f 0.3 -s -wao -a - -b data_miRNAs/9mirnas/Pre-miRNA.gff3 | awk -F"\t" '$7!="."' | cut -f 4,15 \
#   | sed -r "s/\t.+Name=/\t/" | awk '{OFS="\t";if(a[$2]){print $1,$2a[$2]} else {a[$2]++;print $1,$2}}' >data_miRNAs/9mirnas/NameChange.txt

# awk -F"\t" 'NR==FNR{a[$1]=$2;next}{OFS="\t";if(a[$22]){$22=a[$22]};print $0}' data_miRNAs/9mirnas/NameChange.txt data_miRNAs/9mirnas/tmpTable1.txt > data_miRNAs/9mirnas/tmp_table.txt
# sort -k 9,9 data_miRNAs/9mirnas/tmp_table.txt | awk '{print ">"""$22"""_"""$10"""\n"""$9}'  >data_miRNAs/9mirnas/pre_name_seq.fa

# ssearch36 -m 8 -E 1 data_miRNAs/9mirnas/pre_name_seq.fa data_miRNAs/9mirnas/pre_name_seq.fa > data_miRNAs/9mirnas/2result.txt
# awk -F"\t" '{OFS="\t";split($1,c,"_");split($2,d,"_");a=c[2];b=int($3*$4/100+0.5-$6);if(a-b<=2&&c[1]!=d[1]) print c[1],d[1]}' data_miRNAs/9mirnas/2result.txt \
#  | sort -k 1,1 -k 2,2 | uniq >data_miRNAs/9mirnas/Conserved3.txt

# cut -f 1 data_miRNAs/9mirnas/tmpTable1.txt | awk -F"\t" 'NR==FNR{a[$1]=$1;a[$2]=$2;next}{if(!a[$1]) print $1}' data_miRNAs/9mirnas/Conserved3.txt - >data_miRNAs/9mirnas/non_Conserved3.txt

# grep ">" data_miRNAs/9mirnas/pre_name_seq.fa | sed -r 's/>//;s/_..//' | awk -F"\t" 'NR==FNR{a[$1]=$1;a[$2]=$2;next}{if(!a[$1]) print $1}' data_miRNAs/9mirnas/Conserved3.txt - >data_miRNAs/9mirnas/non_Conserved3.txt


# Rscript software_scripts/06-result-table.R

# awk -F"\t" 'NR==FNR{a[$2]=$1;next}{OFS="\t";if(a[$1]){print a[$1],$1}else{print $0}}' data_miRNAs/9mirnas/NameChange.txt  data_miRNAs/9mirnas/NameChange2.txt  >data_miRNAs/9mirnas/NameChange_final.txt

# python software_scripts/06-result-table.py data_miRNAs/9mirnas/mir_seq_len.txt data_miRNAs/8bias/3-pipe_pre-miRNA.str \
# data_miRNAs/8bias/1-pipeline_result.txt data_miRNAs/7isomiRs/alignment_tpm_result.txt data_miRNAs/8bias/2-bias_structure_result.txt  \
# data_miRNAs/9mirnas/tmpTable2.txt data_miRNAs/9mirnas/ Zea_mays data_miRNAs/9mirnas/NameChange_final.txt;

# awk -F"\t" '$11~/3p|5p/&&$18~/3p|5p/' data_miRNAs/9mirnas/tmpTable2.txt | sort -k 23,23 -k 22,22 | awk '{OFS="\t";gsub("MI","mi",$8);gsub("MI","mi",$15);print $0}' > data_miRNAs/9mirnas/finalTable.txt


# python software_scripts/7_Hairpin_plots.py data_miRNAs/10plot data_miRNAs/9mirnas/finalTable.txt data_miRNAs/7isomiRs/alignment_tpm_result.txt data_miRNAs/6mireapout/posMirnaLocation.txt

pdfunite $(ls -v data_miRNAs/10plot/*.pdf) FigureS1.pdf
