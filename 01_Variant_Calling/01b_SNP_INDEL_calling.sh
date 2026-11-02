#!/bin/bash
# ==============================================================================
# Script: 01b_Illumina_SNP_calling.sh
# Purpose: GATK4 Best Practices Pipeline for Population-Scale SNP/INDEL Calling
# Description: Executes Round 1 calling, Hard Filtering, BQSR, and Round 2 calling.
# Author: Zhenzhen Sun
# Target: High-altitude Pig Adaptation Project
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. Environment & Paths
# ------------------------------------------------------------------------------
source activate /storage/public/home/rengang/miniconda3/envs/gatk

GATK_JAR="/storage/public/home/rengang/001liu/project/epi/data/szz/soft/gatk/build/gatk-4.6.0.0-8-g64dc4b3-SNAPSHOT/gatk-package-4.6.0.0-8-g64dc4b3-SNAPSHOT-local.jar"
REF="/storage/public/home/rengang/001liu/project/epi/data/szz/ref_seq/GCF_000003025.6_Sscrofa11.1_genomic.fna"
TMPDIR="/storage/public/home/rengang/001liu/project/epi/data/szz/TMP"
BAM_DIR="/storage/public/home/rengang/001liu/project/epi/data/szz/BAM_FILES"
INTERVAL_BED="/storage/public/home/rengang/001liu/project/epi/data/szz/test/undo.bed"
THREADS=36

# Create a temporary directory
mkdir -p ${TMPDIR}

log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

# ------------------------------------------------------------------------------
# 2. Merge vcf -> Separate SNPs/INDELS -> Calculate DP -> Hard filtering
# ------------------------------------------------------------------------------
process_and_filter_vcf() {
    local ROUND_NAME=$1
    log "Starting VCF processing and filtering for: ${ROUND_NAME}"
    
    mkdir -p vcf-files_${ROUND_NAME} idx-files_${ROUND_NAME}
    
    mv *.vcf vcf-files_${ROUND_NAME}/ 2>/dev/null || true
    mv *.idx idx-files_${ROUND_NAME}/ 2>/dev/null || true
    
    cd vcf-files_${ROUND_NAME}/
    
    log "Merging scattered VCFs..."
    onevcf=$(ls *.vcf | head -n 1)
    grep "^#" ${onevcf} > head.txt
    for f in *.vcf; do grep -v "^#" $f >> body.txt; done
    cat head.txt body.txt > ../merged_${ROUND_NAME}.vcf
    rm head.txt body.txt
    cd ..

    log "Separating SNPs and INDELs..."
    vcftools --vcf merged_${ROUND_NAME}.vcf --remove-indels --recode --recode-INFO-all --out ${ROUND_NAME}_snps-only
    vcftools --vcf merged_${ROUND_NAME}.vcf --keep-only-indels --recode --recode-INFO-all --out ${ROUND_NAME}_indels-only
    
    log "Sorting VCFs..."
    cat ${ROUND_NAME}_snps-only.recode.vcf | vcf-sort -t $TMPDIR -p ${THREADS} -c > ${ROUND_NAME}_sorted-snps-only.vcf
    cat ${ROUND_NAME}_indels-only.recode.vcf | vcf-sort -t $TMPDIR -p ${THREADS} -c > ${ROUND_NAME}_sorted-indels-only.vcf

    log "Calculating dynamic MAX_DEPTH for filtering..."
    grep -v "^#" ${ROUND_NAME}_snps-only.recode.vcf | cut -f 8 | grep -oe ";DP=.*" | cut -f 2 -d ";" | cut -f 2 -d "=" > dp_snps.txt
    MAXDEPTH=$(cat dp_snps.txt | grep -v "^N" | sort -n | awk '{print ($(NF-1)+$NF)*5}')
    log "Calculated MAX_DEPTH for SNPs: ${MAXDEPTH}"

    log "Applying GATK Hard Filtering on SNPs..."
    java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar ${GATK_JAR} VariantFiltration \
        --reference ${REF} \
        --variant ${ROUND_NAME}_sorted-snps-only.vcf \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 45.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > ${MAXDEPTH}" \
        --filter-name "FAIL" \
        --output ${ROUND_NAME}_sorted-snps_filtered.vcf

    log "Applying GATK Hard Filtering on INDELs..."
    java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar ${GATK_JAR} VariantFiltration \
        --reference ${REF} \
        --variant ${ROUND_NAME}_sorted-indels-only.vcf \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "FAIL" \
        --output ${ROUND_NAME}_sorted-indels_filtered.vcf

    log "Extracting PASS-only SNPs for downstream processes..."
    grep "^#" ${ROUND_NAME}_sorted-snps_filtered.vcf > head.vcf
    grep -v "^#" ${ROUND_NAME}_sorted-snps_filtered.vcf | awk '$7=="PASS"' > pass.vcf
    cat head.vcf pass.vcf > ${ROUND_NAME}_sorted-snps_filtered-pass-only.vcf
    rm head.vcf pass.vcf
    
    log "Completed processing for ${ROUND_NAME}!"
}

# ------------------------------------------------------------------------------
# 3. ROUND 1: Initial HaplotypeCaller (Generation of Raw VCF)
# ------------------------------------------------------------------------------
log "================ STARTING ROUND 1 CALLING ================"
BAM_INPUTS_R1=$(find ${BAM_DIR} -name "*_final_new.bam" -printf "-I %p ")

while read line; do
    out=$(echo $line | sed 's/:/_/g')
    log "Calling variants on interval: $line"
    java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
        -R ${REF} ${BAM_INPUTS_R1} -L $line -O ${out}.vcf
done < ${INTERVAL_BED}

process_and_filter_vcf "first-round"
KNOWN_SITES_VCF="first-round_sorted-snps_filtered-pass-only.vcf"

# ------------------------------------------------------------------------------
# 4. BQSR: Base Quality Score Recalibration
# ------------------------------------------------------------------------------
log "================ STARTING BQSR ================"
for IBAM in $(find ${BAM_DIR} -name "*_final_new.bam"); do
    OBAM=$(basename ${IBAM%.*})
    log "Recalibrating: ${OBAM}"
    
    java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar ${GATK_JAR} BaseRecalibrator \
        --reference ${REF} --input ${IBAM} --known-sites ${KNOWN_SITES_VCF} --output ${BAM_DIR}/${OBAM}_bef-R1.table
        
    java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar ${GATK_JAR} ApplyBQSR \
        --reference ${REF} --input ${IBAM} --output ${BAM_DIR}/${OBAM}_recal.bam --bqsr-recal-file ${BAM_DIR}/${OBAM}_bef-R1.table
done

# ------------------------------------------------------------------------------
# 5. ROUND 2: Final HaplotypeCaller on Recalibrated BAMs
# ------------------------------------------------------------------------------
log "================ STARTING ROUND 2 CALLING ================"
BAM_INPUTS_R2=$(find ${BAM_DIR} -name "*_recal.bam" -printf "-I %p ")

while read line; do
    out=$(echo $line | sed 's/:/_/g')
    log "Calling variants on interval: $line (Recalibrated)"
    java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar ${GATK_JAR} HaplotypeCaller \
        -R ${REF} ${BAM_INPUTS_R2} -L $line -O ${out}_round2.vcf
done < ${INTERVAL_BED}

process_and_filter_vcf "second-round"

log "================ ALL DONE SUCCESSFULLY ================"
