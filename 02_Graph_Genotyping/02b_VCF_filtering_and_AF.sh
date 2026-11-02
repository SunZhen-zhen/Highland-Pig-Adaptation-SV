#!/bin/bash
# ==============================================================================
# Script: 02b_VCF_filtering_and_AF.sh
# Purpose: Rigorous Quality Control and Allele Frequency Extraction
# Description: Filters out low-confidence graph calls and calculates cohort-specific AF.
# Author: Zhenzhen Sun
# Target: High-altitude Pig Adaptation Project
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. Environment & Paths
# ------------------------------------------------------------------------------
source activate /storage/public/home/rengang/sunzhenzhen/miniconda3_new/envs/NGS

# Paragraph The original VCF file generated from running (bgzf format)
RAW_VCF="population_sv_merged_bgzf.vcf.gz"
TAGGED_VCF="population_sv_tagged.vcf.gz"
FINAL_VCF="population_sv_clean.vcf.gz"
AF_OUTPUT="population_frequencies.txt"

# File containing group classification information, used for extracting frequencies by group
# The format is usually: sample_name population_name
SAMPLE_SPECIES_FILE="sample_species.txt"

log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

# ------------------------------------------------------------------------------
# 2. Calculating Tags and Metrics
# ------------------------------------------------------------------------------
log "================ CALCULATING POPULATION TAGS ================"
# Use the "fill-tags" function to fill in the missing INFO fields (such as basic information like AF, AC, AN, etc.)
bcftools +fill-tags ${RAW_VCF} -O z -o ${TAGGED_VCF} -- -t all
bcftools index -t ${TAGGED_VCF}

# ------------------------------------------------------------------------------
# 3. Stringent Call-rate & Quality Filtering (As described in Methods)
# ------------------------------------------------------------------------------
log "================ EXECUTING STRINGENT FILTERING ================"
# Remove the loci that fail quality control:
# 1. F_MISSING < 0.2 (The rate of missing data within the population should not exceed 20%)
# 2. QU >= 10 (The genotyping quality provided by Paragraph must reach an extremely high level of confidence)
bcftools filter \
  -i 'F_MISSING < 0.2 && QU >= 10' \
  ${TAGGED_VCF} \
  -O z \
  -o ${FINAL_VCF}

bcftools index -t ${FINAL_VCF}

RETAINED_SVS=$(bcftools stats ${FINAL_VCF} | grep "number of records:" | awk -F ':' '{print $2}' | xargs)
log "Filtration complete! High-confidence non-redundant SVs retained: ${RETAINED_SVS}"

# ------------------------------------------------------------------------------
# 4. Extracting Cohort-Specific Allele Frequencies
# ------------------------------------------------------------------------------
log "================ EXTRACTING ALLELE FREQUENCIES ================"
# It is necessary to use the -S option to specify the external sample file so that bcftools can calculate the AF (allele frequency) for different groups separately.
bcftools +fill-tags ${FINAL_VCF} -- -S ${SAMPLE_SPECIES_FILE} -t AF | \
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF_Huzhu_Bamei\t%INFO/AF_Huangzhong_Bamei\t%INFO/AF_Hanjiang\t%INFO/AF_Changbaishan_WildPig\t%INFO/AF_Tibetan\n' \
> ${AF_OUTPUT}

log "Allele frequencies successfully extracted to ${AF_OUTPUT}."
log "Pipeline 02 execution completed!"
