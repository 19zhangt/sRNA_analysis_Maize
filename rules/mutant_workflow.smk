
# Time: 2019-10
# sRNA-seq pipeline with snakemake

shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

OUTPUTLIST = []

OUTPUTLIST.extend(expand("3mapdata/{sample}.bw", sample=config["samples"]))
OUTPUTLIST.extend(expand("4readlength/{sample}_all.txt", sample=config["samples"]))
OUTPUTLIST.extend(expand("4rpmdata/{sample}.txt", sample=config["samples"]))
OUTPUTLIST.extend(expand("3mapdata/{sample}.bw", sample=config["samples"]))
# OUTPUTLIST.extend(["4readlength/"])

localrules: all
rule all:
    input: OUTPUTLIST

rule raw_data_gz:
    input: "0rawdata/{sample}.fastq.gz"
    output: "0unzip/{sample}.fastq"
    threads: 5
    shell: "gunzip -c {input} > {output}"

rule read_clean:
    input: "0unzip/{sample}.fastq"
    output: "1cleandata/{sample}_clean.fasta"
    log: "log/0adapter_phred.txt"
    threads: 5
    shell:
        """
        adapter=`grep {wildcards.sample} Adapter.txt | cut -f 2`
        ../software_scripts/fastx/fastx_clipper -n -v -Q33 -l 10 -a $adapter -i {input} -o {input}.tmp1
        ../software_scripts/fastx/fastq_quality_filter -q 20 -p 80 -Q33 -i {input}.tmp1 -o {input}.tmp2
        ../software_scripts/fastx/fastq_to_fasta -r -n -v -Q33 -i {input}.tmp2 -o {output}
        rm {input}.tmp1 {input}.tmp2
        """


rule length_filter:
    input: "1cleandata/{sample}_clean.fasta"
    output: "1cleandata/{sample}.fasta"
    threads: 5
    shell:
        """
        set -x
        perl ../software_scripts/perl_length_cutoff.pl -i {input} -min 18 -max 26 -name {output}
        """

rule collapse_fa:
    input: "1cleandata/{sample}.fasta"
    output: "2collapsedata/{sample}.fasta"
    threads: 5
    shell:
        """
        ../software_scripts/fastx/fastx_collapser -v -i {input} -o {output}
        """


rule read_map:
    input: "2collapsedata/{sample}.fasta"
    output:
        bw = "3mapdata/{sample}.bw",
        gm = "3mapdata/{sample}_genome.fasta",
        clm = "3mapdata/{sample}_rmtrsno.fasta"
    threads: 5
    log: "log/00Table_Summary_of_sRNA-seq_data_in_differnet_samples.txt"
    shell:
        """
        bowtie -p 5 -v 1 -f -t -m 20 -S --al {output.gm} ../genome/Genome/Genome {input} >3mapdata/{wildcards.sample}.sam
        samtools view -b -S 3mapdata/{wildcards.sample}.sam > 3mapdata/{wildcards.sample}.bam
        samtools sort -@ 6  3mapdata/{wildcards.sample}.bam -o 3mapdata/{wildcards.sample}_sort.bam
        samtools index -@ 6 3mapdata/{wildcards.sample}_sort.bam
        bamCoverage -b 3mapdata/{wildcards.sample}_sort.bam -o {output.bw} --normalizeUsing None -bs 5 -p 5
        rm 3mapdata/{wildcards.sample}.*am
        bowtie -p 5 -v 0 -f -t --norc -S --un {output.clm} ../genome/trsnsnoRNAs/trsnsnoRNAs {output.gm} 1>/dev/null 2>&1
        echo -e "ID\tTotal Reads\tClean Reads\tClean Ratio\tGenome Match\tRemove t/r/sn/snoRNA" >{log}
        """

rule summary_read:
    input:
        fq =  "0unzip/{sample}.fastq",
        clfa = "1cleandata/{sample}.fasta",
        ccfa = "2collapsedata/{sample}.fasta",
        gm = "3mapdata/{sample}_genome.fasta",
        clm = "3mapdata/{sample}_rmtrsno.fasta"
    output:
        all = "4readlength/{sample}_all.txt",
        unique = "4readlength/{sample}_unique.txt"
    log: "log/00Table_Summary_of_sRNA-seq_data_in_differnet_samples.txt"
    shell:
        """
        readnum=$((`cat {input.fq} | wc -l`/4))
        clean=`grep ">" {input.clfa} | wc -l`
        ratioval=`awk 'BEGIN{{printf ('$clean'/'$readnum')*100}}'`
        gg=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input.gm}`
        rmtr=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input.clm}`
        echo -e "{wildcards.sample}\t${{readnum}}\t${{clean}}\t${{ratioval}}\t${{gg}}\t${{rmtr}}" >>{log}
        awk -v sample={wildcards.sample} '$1!~/>/{{a[length($1)]++}}END{{OFS="\t";for(i in a)print "All",sample,i,a[i]}}' {input.clfa} >{output.all}
        awk -v sample={wildcards.sample} '$1!~/>/{{a[length($1)]++}}END{{OFS="\t";for(i in a)print "All",sample,i,a[i]}}' {input.ccfa} >{output.unique}
        """

rule read_rpm:
    input: "3mapdata/{sample}_rmtrsno.fasta"
    output: "4rpmdata/{sample}.txt"
    threads: 5
    shell:
        """
        readmap=`awk '$1~/>/{{split($1,a,"-");b+=a[2]}}END{{print b}}' {input}`
        awk -v RS=">" -v SRR={wildcards.sample} 'NR>1{{split($1, a,"-");RPM=a[2]/'$readmap'*1000000;if(a[2]>5){{print SRR"\t"$2"\t"RPM}}}}' {input} >{output}
        """

rule merge_sample:
    output:
        read = "data_mutant/5largematrix/unique_read_RPM.txt",
        mir = "data_mutant/5largematrix/mireap_20_24.fa",
        mm = "data_mutant/5largematrix/reads_18_26.fa"
    shell:"Rscript software_scripts/sample_merge.R data_mutant/4rpmdata/ {output.read} 100 10 20 24 18 26 {output.mir} {output.mm}"
