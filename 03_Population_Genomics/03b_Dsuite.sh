#!/bin/bash
# ==============================================================================
# Script: 03b_Dsuite.sh
# Purpose: Formal Introgression Testing (ABBA-BABA and f4-ratio)
# Description: Uses unpruned whole-genome SNPs to detect historical gene flow.
# Author: Zhenzhen Sun
# ==============================================================================

set -euo pipefail
source activate /storage/public/home/rengang/miniconda3/envs/raxml-ng

INPUT_VCF="filtered_snps.recode_filtered.vcf"
# This file contains two columns: Sample_ID and Population_Name. (The Changbaishan wild boar must be named as Outgroup)
POP_MAP="sets_dsuite.txt" 
# Mandatory tree topology file to ensure that Dsuite is aware of species relationships
TREE="tree_pre.txt"

echo "[INFO] Running Dsuite Dtrios for ABBA-BABA D-statistics..."
# -t specifies the tree topology, -c specifies the calculation of f4-ratio
Dsuite Dtrios -t ${TREE} -c ${INPUT_VCF} ${POP_MAP}

echo "[INFO] Summarizing D-statistics Z-scores block-jackknifing..."
Dsuite Fbranch ${TREE} sets_dsuite_BBAA.txt > fbranch_results.txt

echo "[INFO] 03b_Dsuite completed successfully. Please check Bamei_Introgression_BBAA.txt for |Z|<3 thresholds."
