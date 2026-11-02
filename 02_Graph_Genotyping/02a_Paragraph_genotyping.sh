#!/bin/bash
# ==============================================================================
# Script: 02a_Paragraph_genotyping.sh
# Purpose: Population-Scale Graph-Based SV Genotyping using Paragraph
# Description: Maps short-read NGS data to the HiFi-derived SV graph backbone.
# Author: Zhenzhen Sun
# Target: High-altitude Pig Adaptation Project
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. Environment & Global Parameters
# ------------------------------------------------------------------------------
source activate /storage/public/home/rengang/sunzhenzhen/miniconda3_new/envs/NGS 
REF="/storage/public/home/rengang/001liu/project/epi/data/szz/ref_seq/GCF_000003025.6_Sscrofa11.1_genomic.fna"

# step 01a backbone file
GRAPH_BACKBONE="pig_sv_graph_backbone.vcf.gz"

# A list file containing 150 sample paths, sample names and BAM paths
# The format is usually: id  path (for example: Bamei_01 /path/to/Bamei_01.bam)
MANIFEST="manifest_accurate.txt"

OUT_DIR="paragraph_output"
MAX_DEPTH=20
THREADS=64

log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

# ------------------------------------------------------------------------------
# 2. Graph Genotyping using Paragraph (multigrmpy.py)
# ------------------------------------------------------------------------------
log "================ STARTING PARAGRAPH GENOTYPING ================"
log "Using SV backbone: ${GRAPH_BACKBONE}"
log "Manifest file: ${MANIFEST}"

# Run Paragraph (multigrmpy.py is the core module for handling groups)
# -M: Maximum read depth parameter to prevent memory exhaustion in highly repetitive centromeric arrays.
multigrmpy.py \
  -i ${GRAPH_BACKBONE} \
  -m ${MANIFEST} \
  -r ${REF} \
  -o ${OUT_DIR} \
  -M ${MAX_DEPTH} \
  --threads ${THREADS}

log "Paragraph genotyping completed successfully. Output stored in ${OUT_DIR}."
# Note: After the Paragraph process is completed, the merged original genotype results are usually located at ${OUT_DIR}/genotypes.vcf.gz
