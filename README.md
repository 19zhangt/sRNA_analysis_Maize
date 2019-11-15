## small RNA-seq data analysis in maize (The scripts used in manuscripts)

This repository is being developed for providing an entire workflow for processing small RNA sequencing data in maize. This directory includes key scripts for analyzing miRNAs, isomiRs, phasiRNAs, mutant and degradome data.

## How to start:

### Set up the dependency environment based on conda

```
## Creating an environment
conda env create -f environment.yml

## Exporting an environment
conda env export > environment.yml
```



### Run the pipeline based on snakemake

```
## Drawing the scheme
snakemake --snakefile sRNA_workflow.smk -np --dag | dot -Tpdf > sRNA_workflow.pdf

## Running the pipeline
snakemake --snakefile sRNA_workflow.smk --cores 20 --latency-wait 120
```

