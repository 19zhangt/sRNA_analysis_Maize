
# Time: 2019-11
# degradome data pipeline with snakemake

rule raw_data_gz:
    input: "data_degradome/0rawdata/{sample}.fastq.gz"
    output: "data_degradome/0unzip/{sample}.fastq"
    threads: 5
    shell: "gunzip -c {input} > {output}"

# rule read_clean:
#     input: "data_degradome/0unzip/{sample}.fastq"
#     output: "data_degradome/1cleandata/{sample}_clean.fasta"
#     threads: 5
#     shell:
#         """
#         adapter=`grep {wildcards.sample} Adapter.txt | cut -f 2`
#         software_scripts/fastx/fastx_clipper -n -v -Q33 -l 10 -a $adapter -i {input} -o {input}.tmp1
#         software_scripts/fastx/fastq_quality_filter -q 20 -p 80 -Q33 -i {input}.tmp1 -o {input}.tmp2
#         software_scripts/fastx/fastq_to_fasta -r -n -v -Q33 -i {input}.tmp2 -o {output}
#         rm {input}.tmp1 {input}.tmp2
#         """
#
# rule length_filter:
#     input: "data_degradome/1cleandata/{sample}_clean.fasta"
#     output: "data_degradome/1cleandata/{sample}.fasta"
#     threads: 5
#     shell:
#         """
#         perl software_scripts/perl_length_cutoff.pl -i {input} -min 20 -max 21 -name {output}
#         """
