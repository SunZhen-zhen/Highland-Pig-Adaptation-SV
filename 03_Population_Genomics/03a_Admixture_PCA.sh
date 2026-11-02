#!/bin/bash
# ==============================================================================
# Script: 03a_Admixture_PCA.sh
# Purpose: Population Structure Analysis (LD Pruning, PCA, and ADMIXTURE)
# Description: Filters SNPs for Hardy-Weinberg, LD, MAF, then models ancestry.
# Author: Zhenzhen Sun
# ==============================================================================

set -euo pipefail
source activate /storage/public/home/rengang/miniconda3/envs/raxml-ng

INPUT_VCF="filtered_snps.recode_filtered.vcf"
THREADS=36

echo "[INFO] Starting comprehensive Quality Control using PLINK..."
plink --vcf ${INPUT_VCF} \
      --maf 0.05 \
      --mind 0.1 \
      --geno 0.1 \
      --hwe 1e-6 \
      --allow-extra-chr \
      --make-bed \
      --out qc_data

echo "[INFO] Performing Linkage Disequilibrium (LD) pruning..."
plink --bfile qc_data \
      --indep-pairwise 50 10 0.2 \
      --allow-extra-chr \
      --out pruned_data

echo "[INFO] Extracting pruned neutral SNPs..."
plink --bfile qc_data \
      --extract pruned_data.prune.in \
      --allow-extra-chr \
      --make-bed \
      --out pruned_snps

echo "[INFO] Running Principal Component Analysis (PCA)..."
plink --bfile pruned_snps \
      --pca \
      --allow-extra-chr \
      --out SNP_PCA_results

echo "[INFO] Running ADMIXTURE cross-validation for K=2 to 6..."
for K in {2..6}; do
    echo "Running K=${K}..."
    admixture --cv pruned_snps.bed ${K} -j${THREADS} | tee log_admixture_${K}.out
done

echo "[INFO] 03a_Admixture_PCA completed successfully."
