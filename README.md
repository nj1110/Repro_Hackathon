# Hackaton project
## _Group 11: Perform analysis of diffential slicing based on RNA-Seq samples from patients with uveal melanoma_ 
Université Paris Saclay - M2 AMI2B



## Project Description

Objective of this pipeline is to reproduce the analyses led in 2 publications (Harbour et al 2013, Furney et al 2017) about the effect of SF3B1 mutation on alternative slicing of genes coding for proteins potentially involved in uveal melanoma. 

> The whole process consists of donwloading patients RNA-seq SRR as fastq files, 
> ensuring quality conrtrol, ensuring aligned reads and counting, 
> visualizing and assessing alternative slicing,
> using different tools containerized thanks to Docker 
(sra-tools, FastQC, rna-star, subread featureCounts, R and DESeq2 package) 
> within Snakemake workflows

## Organization
- a subfolder with dockerfiles
- three snakefiles with workflows:
    - one worklfow to download data and create index
    - one workflow for mapping and counts  
    - one worflow for XXX
- a README file  
- a RUN file

## Installation and execution

As a prerequisite you need conda, Docker, Singularity installed.
We recommend to use VM with 64 Go RAM to run the whole pipeline.

Please make sure that you have activated the conda environment:

```sh
conda activate snakemake
```

To run the whole process please follow the following instructions below:
```sh
snakemake -s wkf_download-data_8SRR_index.txt --cores 8 --use-singularity
snakemake -s wkf_map-analyze_8SRR.txt --cores 2 --use-singularity
```

or 
```sh
bash run.sh
```





