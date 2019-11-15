
# Time: 2019-10
# sRNA-seq pipeline with snakemake

shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")



# include : "rules/miRNAs_workflow.smk"


configfile: "config.yaml"

include : "rules/phasiRNAs_workflow.smk"

rule all:
    input:
        expand("data_phasiRNA_v2/1phasiRNAs/{sample}.log", sample=config["extract_samples"])

# rule all:
#     input:
#         expand("data_miRNAs/4countdata/{sample}.txt", sample=config["samples"]),
#         expand("data_miRNAs/5largematrix_count/unique_read_count.txt")

# include : "rules/degradome_workflow.smk"

# rule degradome_out:
#     input:
#         expand("data_degradome/0unzip/{sample}.fastq", sample=config["degradome_data"])
