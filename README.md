# Highland-Pig-Adaptation-SV
# Genomic Architecture of High-Altitude Adaptation in Pigs

This repository contains the custom code, bioinformatics pipelines, and analytical scripts used in the study:
Ancient structural variants facilitate high-altitude adaptation in pigs through targeted regulatory and topological restructuring (Currently under review).

## Overview
We deployed a hybrid long-read/short-read graph-pangenome strategy to resolve the population-scale structural variation (SV) landscape across diverse porcine lineages. The scripts provided herein allow for the reproduction of variant discovery, graph-based genotyping, and downstream evolutionary genomics analyses.

## Directory Structure
* `01_Variant_Calling/`: Scripts for PacBio HiFi SV discovery (Sniffles2) and Illumina SNP/INDEL calling (GATK).
* `02_Graph_Genotyping/`: Population-scale SV genotyping using Paragraph and rigorous VCF filtering.
* `03_Population_Genomics/`: Scripts for ADMIXTURE, PCA, D-statistics, and selective sweep ($Z-F_{ST}$, $\pi$) analyses.
* `04_Borzoi_Analysis/`: Custom scripts for mechanism elucidation and multi-omics visualization.

## Software Prerequisites
Major software utilized in this study includes (but is not limited to):
* PacBio pbmm2 (v26.2.0) & Sniffles2 (v2.4)
* Paragraph (v2.4a)
* BWA-MEM (v0.7.17) & GATK (v4.7.0)
* BCFtools (v1.25) & PLINK (v1.9)
* ADMIXTURE (v1.3.0) & TreeMix (v1.13)
* Dsuite (v0.5)


## Data Availability
The raw sequencing data generated in this study are available in the NCBI Sequence Read Archive (SRA) under BioProject PRJNA1550252. 

## Contact
For questions regarding the code or bioinformatic pipelines, please contact: 2025060294@nwafu.edu.cn or open an Issue in this repository.
