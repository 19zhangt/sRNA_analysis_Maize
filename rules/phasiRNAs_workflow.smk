
# Time: 2019-10
# phasiRNA pipeline with snakemake

# rule length_filter:
#     input: "data_miRNAs/1cleandata/{sample}.fasta"
#     output:
#         tofile = "data_phasiRNA_v2/1cleandata/{sample}_21.fasta",
#         tffile = "data_phasiRNA_v2/1cleandata/{sample}_24.fasta"
#     threads: 10
#     shell:
#         """
#         perl software_scripts/perl_length_cutoff.pl -i {input} -min 21 -max 21 -name {output.tofile}
#         perl software_scripts/perl_length_cutoff.pl -i {input} -min 24 -max 24 -name {output.tffile}
#         """

rule phasiRNA:
    input:
        tofile = "data_phasiRNA_v2/1cleandata/{sample}_21.fasta",
        tffile = "data_phasiRNA_v2/1cleandata/{sample}_24.fasta"
    output: "data_phasiRNA_v2/1phasiRNAs/{sample}.log"
    threads: 10
    shell:
        """
        codePath=`pwd`
        export PERL5LIB=$PERL5LIB:$codePath/software_scripts/Perl_modules/share/perl/5.18.2
        python3 software_scripts/phasdetect.py index/Genome_1_10 data_phasiRNA_v2/1phasiRNAs/{wildcards.sample} 21 {input.tofile} 2>{output}
        python3 software_scripts/phasdetect.py index/Genome_1_10 data_phasiRNA_v2/1phasiRNAs/{wildcards.sample} 24 {input.tffile} 2>{output}
        """
