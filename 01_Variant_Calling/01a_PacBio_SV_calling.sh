#!/bin/bash
# ==============================================================================
# Script: 01a_PacBio_SV_calling.sh
# Purpose: PacBio HiFi Alignment and De Novo SV Discovery (Graph Backbone Prep)
# Description: pbmm2 alignment, Sniffles2 joint-calling, and stringent filtering.
# Author: ZhenzhenSun
# Target: High-altitude Pig Adaptation Project
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. Environment & Paths
# ------------------------------------------------------------------------------
source activate /storage/public/home/rengang/miniconda3/envs/sv_calling
REF="/storage/public/home/rengang/001liu/project/epi/data/szz/ref_seq/GCF_000003025.6_Sscrofa11.1_genomic.fna"
RAW_DIR="/storage/public/home/rengang/001liu/project/epi/data/szz/Pig_ear_three_generations_HiFi/raw_reads"
BAM_DIR="/storage/public/home/rengang/001liu/project/epi/data/szz/Pig_ear_three_generations_HiFi/bam_sorted_file"
THREADS=36

mkdir -p ${BAM_DIR}

log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

# ------------------------------------------------------------------------------
# 2. Alignment using pbmm2
# ------------------------------------------------------------------------------
log "================ STARTING HIFI ALIGNMENT ================"
for FASTQ in $(ls ${RAW_DIR}/*.fastq.gz); do
    NAME=$(basename ${FASTQ%.fastq.gz})
    log "Aligning sample: ${NAME}"
    
    pbmm2 align ${REF} ${FASTQ} ${BAM_DIR}/${NAME}_sorted.bam \
        --preset CCS \
        --sort \
        -j ${THREADS}
        
    samtools index -@ ${THREADS} ${BAM_DIR}/${NAME}_sorted.bam
done

# ------------------------------------------------------------------------------
# 3. SV Discovery: Sniffles2 Single-Sample Feature Extraction (.snf)
# ------------------------------------------------------------------------------
log "================ STARTING SNIFFLES2 FEATURE EXTRACTION ================"
mkdir -p ${BAM_DIR}/snf_files
rm -f snf_files_list.tsv

for BAM in $(ls ${BAM_DIR}/*_sorted.bam); do
    NAME=$(basename ${BAM%_sorted.bam})
    log "Extracting SV features for: ${NAME}"
    
    sniffles --input ${BAM} \
             --snf ${BAM_DIR}/snf_files/${NAME}.snf \
             --threads ${THREADS}
             
    echo "${BAM_DIR}/snf_files/${NAME}.snf" >> snf_files_list.tsv
done

# ------------------------------------------------------------------------------
# 4. SV Joint Calling & Merging (The 'Mosaic' Discovery Panel)
# ------------------------------------------------------------------------------
log "================ STARTING SNIFFLES2 JOINT CALLING ================"
sniffles --input snf_files_list.tsv \
         --mosaic \
         --threads ${THREADS} \
         --vcf pig_sv_merged_raw.vcf

# ------------------------------------------------------------------------------
# 5. Stringent SV Filtering & Graph Backbone Preparation
# ------------------------------------------------------------------------------
log "================ PREPARING PRISTINE GRAPH BACKBONE ================"

# Step 5.1: Basic Filtering (Keep only PASS, with a length of at least 50 bp, and exclude BND translocations that are not supported by Paragraph)
bcftools view -f PASS pig_sv_merged_raw.vcf | \
    bcftools filter -i 'ABS(SVLEN)>=50 && SVTYPE!="BND"' \
    -O z -o pig_sv_filtered.vcf.gz

bcftools index -t pig_sv_filtered.vcf.gz

# Step 5.2: Remove unnecessary FORMAT fields (required for the standardized input of Graph-genotyper Paragraph)
# Eliminate the individual's genotype details and only provide the pure variation coordinates as the backbone
bcftools annotate -x FORMAT \
    pig_sv_filtered.vcf.gz \
    -O z -o pig_sv_graph_backbone.vcf.gz

bcftools index -t pig_sv_graph_backbone.vcf.gz

log "Completed! The structural backbone (pig_sv_graph_backbone.vcf.gz) is ready for Paragraph genotyping."
