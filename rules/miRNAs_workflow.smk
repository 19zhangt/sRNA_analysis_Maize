
# Time: 2019-10
# sRNA-seq pipeline with snakemake

# OUTPUTLIST = []
# # OUTPUTLIST.extend(expand("data_miRNAs/4readlength/{sample}_all.txt", sample=config["samples"]))
# # OUTPUTLIST.extend(expand("data_miRNAs/4readlength/{sample}_unique.txt", sample=config["samples"]))
# OUTPUTLIST.extend(["data_miRNAs/6mireapout/mireap.txt"])
# # OUTPUTLIST.extend(["data_miRNAs/6mireapout/mireap.txt"])
# # OUTPUTLIST.extend(expand("data_miRNAs/10phasiRNAs/{extract_samples}.log", extract_samples=config["extract_samples"]))
#
# # OUTPUTLIST.extend(expand("data_miRNAs/3mapdata_bw/{sample}.bw", sample=config["samples"]))
# # OUTPUTLIST.extend(expand("data_miRNAs/4rpmdata/{sample}.txt", sample=config["samples"]))
#
# localrules: all
# rule all:
#     input: OUTPUTLIST

# # TGACAGAAGAGAGTGAGCAC
# # PATH=$PATH:`pwd`/software_scripts/sratoolkit/bin:`pwd`/software_scripts/enaBrowserTools/python
# #
# # cat wheat_98 | while read i
# # do
# # enaDataGet -f fastq ${i}
# # if [ ! -f ${i}.fastq.gz ];then
# # 	enaDataGet -f sra ${i}
# # if [ ! -f ${i} ];then
# # 	wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0:6}/${i}/${i}.sra
# # fi
# # fi
# # done
#
# # rule raw_data_sra:
# #     input: "data_miRNAs/0rawdata/{sample}.sra"
# #     output: "data_miRNAs/0unzip/{sample}.fastq"
# #     threads: 5
# #     shell: "./software_scripts/sratoolkit/bin/fastq-dump --split-3 -O data/0unzip/ {input}"
# #
# # rule raw_data_gz:
# #     input: "data_miRNAs/0rawdata/{sample}.fastq.gz"
# #     output: "data_miRNAs/0unzip/{sample}.fastq"
# #     threads: 5
# #     shell: "gunzip -c {input} > {output}"
#
# # rule read_clean:
# #     input: "data_miRNAs/0unzip/{sample}.fastq"
# #     output: "data_miRNAs/1cleandata/{sample}_clean.fasta"
# #     log: "log/0adapter_phred.txt"
# #     threads: 5
# #     shell: "sh software_scripts/read_clean.sh {input} {output} {wildcards.sample} {log}"
#
#
# rule length_filter:
#     input: "data_miRNAs/1cleandata/{sample}_clean.fasta"
#     output: "data_miRNAs/1cleandata/{sample}.fasta"
#     threads: 5
#     shell:
#         """
#         set -x
#         perl software_scripts/perl_length_cutoff.pl -i {input} -min 18 -max 26 -name {output}
#         """
#
# rule collapse_fa:
#     input: "data_miRNAs/1cleandata/{sample}.fasta"
#     output: "data_miRNAs/2collapsedata/{sample}.fasta"wildcards
#     threads: 5
#     shell:
#         """
#         software_scripts/fastx/fastx_collapser -v -i {input} -o {output}
#         """
#
#
# rule read_map:
#     input: "data_miRNAs/2collapsedata/{sample}.fasta"
#     output:
#         gm = "data_miRNAs/3mapdata/{sample}_genome.fasta",
#         clm = "data_miRNAs/3mapdata/{sample}_rmtrsno.fasta",
#         bam = "data_miRNAs/3mapdata/{sample}.bam",
#         bw = "data_miRNAs/3mapdata/{sample}.bw"
#     threads: 5
#     log: "log/00Table_Summary_of_sRNA-seq_data_in_differnet_samples.txt"
#     shell:
#         """
#         bowtie -p 5 -v 1 -f -t -m 20 -S --al {output.gm} genome/Genome/Genome {input} >{output.bam}
#         bamCoverage -b {output.bam} -o {output.bw} --normalizeUsing None -bs 5 -p 5
#         bowtie -p 5 -v 0 -f -t --norc -S --un {output.clm} genome/trsnsnoRNAs/trsnsnoRNAs {output.gm} 1>/dev/null 2>&1
#         echo -e "ID\tTotal Reads\tClean Reads\tClean Ratio\tGenome Match\tRemove t/r/sn/snoRNA" >{log}
#         """

# rule read_bw:
#     input: "data_miRNAs/2collapsedata/{sample}.fasta"
#     output:
#         gm = "data_miRNAs/3mapdata_bw/{sample}_genome.fasta",
#         bw = "data_miRNAs/3mapdata_bw/{sample}.bw"
#     threads: 5
#     shell:
#         """
#         bowtie -p 5 -v 1 -f -t -m 20 -S --al {output.gm} genome/Genome/Genome {input} >data/3mapdata_bw/{wildcards.sample}.sam
#         samtools view -b -S data/3mapdata_bw/{wildcards.sample}.sam > data/3mapdata_bw/{wildcards.sample}.bam
#         samtools sort -@ 6  data/3mapdata_bw/{wildcards.sample}.bam -o data/3mapdata_bw/{wildcards.sample}_sort.bam
#         samtools index -@ 6 data/3mapdata_bw/{wildcards.sample}_sort.bam
#         bamCoverage -b data/3mapdata_bw/{wildcards.sample}_sort.bam -o {output.bw} --normalizeUsing None -bs 5 -p 5
#         rm data/3mapdata_bw/{wildcards.sample}.*am
#         """

# rule summary_read:
#     input:
#         fq =  "data_miRNAs/0unzip/{sample}.fastq",
#         clfa = "data_miRNAs/1cleandata/{sample}.fasta",
#         ccfa = "data_miRNAs/2collapsedata/{sample}.fasta",
#         gm = "data_miRNAs/3mapdata/{sample}_genome.fasta",
#         clm = "data_miRNAs/3mapdata/{sample}_rmtrsno.fasta"
#     output:
#         all = "data_miRNAs/4readlength/{sample}_all.txt",wildcards
#         unique = "data_miRNAs/4readlength/{sample}_unique.txt"
#     log: "log/00Table_Summary_of_sRNA-seq_data_in_differnet_samples.txt"
#     shell:
#         """
#         readnum=$((`cat {input.fq} | wc -l`/4))
#         clean=`grep ">" {input.clfa} | wc -l`
#         ratioval=`awk 'BEGIN{{printf ('$clean'/'$readnum')*100}}'`
#         gg=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input.gm}`
#         rmtr=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input.clm}`
#         echo -e "{wildcards.sample}\t${{readnum}}\t${{clean}}\t${{ratioval}}\t${{gg}}\t${{rmtr}}" >>{log}
#         awk -v sample={wildcards.sample} '$1!~/>/{{a[length($1)]++}}END{{OFS="\t";for(i in a)print "All",sample,i,a[i]}}' {input.clfa} >{output.all}
#         awk -v sample={wildcards.sample} '$1!~/>/{{a[length($1)]++}}END{{OFS="\t";for(i in a)print "All",sample,i,a[i]}}' {input.ccfa} >{output.unique}
#         """
#
# rule summit_length:
#     input: "data_miRNAs/4readlength/"
#     output: "data_miRNAs/4readlength/"
#     log: "log/11Samples_read_summit.txt"
#     shell:
#         """
#         ls {input} | cut -d"_" -f 1 | sort -u | while read line
#         do
#         alln=`awk '$3>17&&$3<27{{print $0}}' ${line}_all.txt | sort -k 4,4n | tail -n 1 | cut -f 3`
#         uniq=`awk '$3>17&&$3<27{{print $0}}' ${line}_unique.txt | sort -k 4,4n | tail -n 1 | cut -f 3`
#         echo -e "{line}\t${alln}\t${uniq}" >>{log}
#         done
#         """
#
#
# rule read_rpm:
#     input: "data_miRNAs/3mapdata/{sample}_rmtrsno.fasta"
#     output: "data_miRNAs/4rpmdata/{sample}.txt"
#     threads: 5
#     shell:
#         """
#         readmap=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input}`
#         awk -v RS=">" -v SRR={wildcards.sample} 'NR>1{{split($1, a,"-");RPM=a[2]/'$readmap'*1000000;if(a[2]>5){{print SRR"\t"$2"\t"RPM}}}}' {input} >{output}
#         """

# rule merge_sample:
#     output:
#         read = "data_miRNAs/5largematrix/unique_read_RPM.txt",
#         mir = "data_miRNAs/5largematrix/mireap_20_24.fa",
#         mm = "data_miRNAs/5largematrix/reads_18_26.fa"
#     shell:"Rscript software_scripts/sample_merge.R data/4rpmdata/ {output.read} 100 10 20 24 18 26 {output.mir} {output.mm}"

rule read_count:
    input: "data_miRNAs/3mapdata/{sample}_rmtrsno.fasta"
    output: "data_miRNAs/4countdata/{sample}.txt"
    threads: 5
    shell:
        """
        readmap=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input}`
        awk -v RS=">" -v SRR={wildcards.sample} 'NR>1{{split($1, a,"-");if(a[2]>5){{print SRR"\t"$2"\t"a[2]}}}}' {input} >{output}
        """

rule merge_sample:
    output:
        read = "data_miRNAs/5largematrix_count/unique_read_count.txt",
        mir = "data_miRNAs/5largematrix_count/mireap_20_24.fa",
        mm = "data_miRNAs/5largematrix_count/reads_18_26.fa"
    shell:"Rscript software_scripts/sample_merge.R data_miRNAs/4countdata/ {output.read} 100 10 20 24 18 26 {output.mir} {output.mm}"


# rule mirep_prepare:
#     input: "data_miRNAs/5largematrix/mireap_20_24.fa"
#     output:
#         text="data_miRNAs/5largematrix/mireap_20_24.txt",
#         sam="data_miRNAs/5largematrix/mireap_20_24.sam"
#     shell:
#         """
#         bowtie -p 10 -f -v 0 -a genome/Genome/Genome {input} >{output.sam}
#         awk -F"\t" '{{OFS="\t";split($1,a," ");print a[1],$3,$4+1,$4+length($5),$2}}' {output.sam} >{output.text}
#         """

# rule mirep_run:
#     input:
#         fasta = "data_miRNAs/5largematrix/mireap_20_24.fa",
#         map = "data_miRNAs/5largematrix/mireap_20_24.txt"
#     output: "data_miRNAs/6mireapout/mireap.txt"
#     shell:
#         """
#         codePath=`pwd`
#         export PERL5LIB=$PERL5LIB:$codePath/software_scripts/mireap_0.2/lib/:$codePath/software_scripts/mireap_0.2/ViennaRNA-2.4.9/interfaces/Perl/
#         export PATH="$PATH:`pwd`/software_scripts/mireap_0.2/ViennaRNA-2.4.9/src/bin"
#         perl software_scripts/mireap_0.2/bin/modified_mireap_loc.pl -i {input.fasta} -m {input.map} -r genome/Genome/Genome.fa -o data/6mireapout/ -A 18 -B 26 -a 18 -b 26 -d 220 -f 20 -p 16
#         python software_scripts/02-MIREAP2isomiRs.py data/6mireapout
#         bash software_scripts/after_mireap.sh
#         """
#

# bedtools intersect -f 0.3 -wao -a genome/miRBasemiRNAs/zma_m.gff3 -b data/6mireapout/mireap-xxx_uniq.gff \
# | grep "miRNA_primary_transcript" | awk -F"\t" '$10=="."' | cut -f 9 | grep -oP "zma.+" \
# | awk -F"\t" 'NR==FNR{a[$1]=$1;next}{if(a[$1]){gsub(/U/,"T",$2);print ">"$1"\n"$2}}' \
# - genome/miRBasemiRNAs/00_miRNAs_sequence.txt >>data/6mireapout/forisoPremirnaSeq.fa
#
# bedtools intersect -f 0.3 -wao -a genome/miRBasemiRNAs/zma_m.gff3 -b data/6mireapout/mireap-xxx_uniq.gff \
# | grep "miRNA_primary_transcript" | awk -F"\t" '$10=="."' | cut -f 9 | grep -oP "zma.+" \
# | awk -F"\t" 'NR==FNR{a[$1]=$1;next}{if(a[$1]){gsub(/U/,"T",$2);print ">"$3"\n"$4"\n>"$5"\n"$6}}' \
# - genome/miRBasemiRNAs/00_miRNAs_sequence.txt >>data/6mireapout/forisoMirnaSeq.fa
#
# bedtools intersect -f 0.3 -wao -a genome/miRBasemiRNAs/zma_m.gff3 -b data/6mireapout/mireap-xxx_uniq.gff \
# | grep "miRNA_primary_transcript" | awk -F"\t" '$10=="."' | cut -f 9 | grep -oP "zma.+" \
# | awk -F"\t" 'NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0}}' \
# - genome/miRBasemiRNAs/01_miRNAs_location.txt >>data/6mireapout/posMirnaLocation.txt
#
# awk -F"\t" 'NR==FNR{a[$1]=$1;next}{if(!a[$1]){gsub(/U/,"T",$2);print ">"$1"\n"$2}}' \
# data/Test.txt genome/miRBasemiRNAs/00_miRNAs_sequence.txt >>data/6mireapout/forisoPremirnaSeq.fa
#
# awk -F"\t" 'NR==FNR{a[$1]=$1;next}{if(!a[$1]){gsub(/U/,"T",$4);print ">"$3"\n"$4"\n>"$5"\n"$6}}' \
# data/Test.txt genome/miRBasemiRNAs/00_miRNAs_sequence.txt >>data/6mireapout/forisoMirnaSeq.fa
#
# awk -F"\t" 'NR==FNR{a[$1]=$0;next}{if(!a[$1]){print $0}}' \
# data/Test.txt genome/miRBasemiRNAs/01_miRNAs_location.txt >>data/6mireapout/posMirnaLocation.txt

# rule phasiRNA:
#     input: "data_miRNAs/1cleandata/{extract_samples}.fasta"
#     output: "data_miRNAs/10phasiRNAs/{extract_samples}.log"
#     threads: 10
#     shell:
        # """
        # codePath=`pwd`
        # export PERL5LIB=$PERL5LIB:$codePath/software_scripts/Perl_modules/share/perl/5.18.2
        # python3 software_scripts/phasdetect.py genome/Genome/Genome data/10phasiRNAs/{wildcards.extract_samples} 21 {input}
        # python3 software_scripts/phasdetect.py genome/Genome/Genome data/10phasiRNAs/{wildcards.extract_samples} 24 {input}
        # """
