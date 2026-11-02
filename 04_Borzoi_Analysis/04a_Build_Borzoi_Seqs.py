#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# Script: 04b_Build_Borzoi_Seqs.py
# Purpose: Construction of 524kb WT and MUT DNA sequences for Deep Learning
# Description: Dynamically parses VCF to extract exact SV breakpoints (INS/DEL), 
#              and reconstructs phased 524kb genomic context for Borzoi inference.
# Author: Zhenzhen Sun
# ==============================================================================

import subprocess
import sys
import os

# ------------------------------------------------------------------------------
# User Configuration
# ------------------------------------------------------------------------------
REF_FASTA    = "GCF_000003025.6_Sscrofa11.1_genomic.fna" # Ensure .fai exists
VCF_FILE     = "population_sv_clean.vcf.gz"              # Bgzipped and indexed VCF
TARGET_SV_ID = "Sniffles2.DEL.2E9M11"                    # Target SV ID (e.g., CNOT2 upstream DEL)
TARGET_LEN   = 524288                                    # Fixed input length required by Borzoi

# ------------------------------------------------------------------------------
# Core Functions
# ------------------------------------------------------------------------------
def run_cmd(cmd):
    """Executes a shell command and returns the output."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {cmd}\n{result.stderr}")
        sys.exit(1)
    return result.stdout.strip()

def get_fasta_seq(chrom, start, end):
    """Extracts reference sequence slice using samtools."""
    cmd = f"samtools faidx {REF_FASTA} {chrom}:{start}-{end}"
    out = run_cmd(cmd)
    lines = out.split('\n')
    seq = "".join(lines[1:]).upper()
    return seq

def get_sv_info(sv_id):
    """Parses exact breakpoint and ALT sequence from the VCF."""
    print(f"[*] Querying VCF for {sv_id} ...")
    cmd = f"bcftools query -i 'ID==\"{sv_id}\"' -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ALT\n' {VCF_FILE}"
    out = run_cmd(cmd)
    
    if not out:
        print(f"[FATAL] SV ID {sv_id} not found in {VCF_FILE}!")
        sys.exit(1)
    
    parts = out.split('\t')
    chrom = parts[0]
    pos = int(parts[1])
    end = int(parts[2]) if parts[2] != "." else pos
    sv_type = parts[3]
    alt_seq = parts[4]
    
    return chrom, pos, end, sv_type, alt_seq

# ------------------------------------------------------------------------------
# Main Logic: Sequence Reconstruction
# ------------------------------------------------------------------------------
def main():
    print(f"============================================================")
    print(f" Constructing Sequences for In Silico Mutagenesis: {TARGET_SV_ID}")
    print(f"============================================================")

    chrom, pos, end, sv_type, alt_seq = get_sv_info(TARGET_SV_ID)
    print(f"    -> Locus: {chrom}:{pos}-{end} ({sv_type})")

    half_target = TARGET_LEN // 2
    wt_out = f"WT_{TARGET_SV_ID}_524k.fa"
    mut_out = f"MUT_{TARGET_SV_ID}_524k.fa"

    if sv_type == "INS":
        ins_len = len(alt_seq)
        if ins_len < 10:
            print("[FATAL] ALT sequence is missing or too short for synthesis.")
            sys.exit(1)
            
        print(f"\n[*] Extracting Wild-Type (WT) Sequence...")
        wt_start = pos - half_target + 1
        wt_end = pos + half_target
        wt_seq = get_fasta_seq(chrom, wt_start, wt_end)
        with open(wt_out, "w") as f:
            f.write(f">{TARGET_SV_ID}_WT_chr{chrom}:{wt_start}-{wt_end}\n{wt_seq}\n")

        print(f"\n[*] Synthesizing Mutant (MUT) Sequence...")
        trim_len = ins_len // 2
        left_flank_len = half_target - trim_len
        right_flank_len = TARGET_LEN - left_flank_len - ins_len

        left_seq = get_fasta_seq(chrom, pos - left_flank_len + 1, pos)
        right_seq = get_fasta_seq(chrom, pos + 1, pos + right_flank_len)

        mut_seq = left_seq + alt_seq + right_seq
        with open(mut_out, "w") as f:
            f.write(f">MUT_{TARGET_SV_ID}_INS{ins_len}\n{mut_seq}\n")

    elif sv_type == "DEL":
        del_len = end - pos + 1
        mid_pos = (pos + end) // 2

        print(f"\n[*] Extracting Wild-Type (WT) Sequence...")
        wt_start = mid_pos - half_target + 1
        wt_end = mid_pos + half_target
        wt_seq = get_fasta_seq(chrom, wt_start, wt_end)
        with open(wt_out, "w") as f:
            f.write(f">{TARGET_SV_ID}_WT_chr{chrom}:{wt_start}-{wt_end}\n{wt_seq}\n")

        print(f"\n[*] Synthesizing Mutant (MUT) Sequence...")
        trim_len = del_len // 2
        left_flank_len = half_target + trim_len
        right_flank_len = TARGET_LEN - left_flank_len

        left_seq = get_fasta_seq(chrom, pos - left_flank_len, pos - 1)
        right_seq = get_fasta_seq(chrom, end + 1, end + right_flank_len)

        mut_seq = left_seq + right_seq
        with open(mut_out, "w") as f:
            f.write(f">MUT_{TARGET_SV_ID}_DEL{del_len}\n{mut_seq}\n")

    else:
        print(f"[ERROR] Unsupported SV type: {sv_type}. Requires INS or DEL.")
        sys.exit(1)

    print(f"\n[+] Success! WT and MUT FASTA files generated.")

if __name__ == "__main__":
    main()
