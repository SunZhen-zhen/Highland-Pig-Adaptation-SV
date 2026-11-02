#!/bin/bash
# ==============================================================================
# Script: 03c_Selective_Sweep.sh
# Purpose: Genome-wide Selection Landscape (Fst, Pi, and Functional Intersection)
# Description: Uses 50-kb sliding windows to identify hard sweeps across lineages.
# Author: Zhenzhen Sun
# ==============================================================================

set -euo pipefail
source activate /storage/public/home/rengang/miniconda3/envs/raxml-ng

INPUT_VCF="use.vcf.gz"
WIN_SIZE=50000
WIN_STEP=25000

echo "[INFO] Calculating Weir and Cockerham F_ST across 50-kb sliding windows..."
vcftools --gzvcf ${INPUT_VCF} \
         --weir-fst-pop HA_list.txt \
         --weir-fst-pop noHA_list.txt \
         --fst-window-size ${WIN_SIZE} \
         --fst-window-step ${WIN_STEP} \
         --out noHA_vs_HA_FST

echo "[INFO] Calculating Nucleotide Diversity (Pi) across 50-kb sliding windows..."
vcftools --gzvcf ${INPUT_VCF} --keep HA_list.txt \
         --window-pi ${WIN_SIZE} --window-pi-step ${WIN_STEP} --out HA_pi

vcftools --gzvcf ${INPUT_VCF} --keep noHA_list.txt \
         --window-pi ${WIN_SIZE} --window-pi-step ${WIN_STEP} --out noHA_pi

echo "[INFO] Integrating functional SVs into selective sweep windows..."
awk 'NR>1 {print $1 "\t" $2-1 "\t" $3}' unique_candidate_windows2.txt > selection_windows.bed

bedtools intersect -a candidate_list_A.vcf -b selection_windows.bed -header > candidate_priority_intersect.vcf

NUM_VARS=$(grep -v "^#" candidate_priority_intersect.vcf | wc -l)
echo "[INFO] Number of prioritized variants residing in hard sweep regions: ${NUM_VARS}"

echo "[INFO] 03c_Selective_Sweep completed successfully."
