#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# Script: 04c_Borzoi_InSilico_Mut.py
# Purpose: Deep Learning (Borzoi) Inference and Publication-ready Visualization
# Description: Executes in silico mutagenesis to predict spatial chromatin 
#              accessibility and RNA-seq output, rendering customized gene tracks.
# Author: Zhenzhen Sun
# ==============================================================================

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

try:
    from baskerville import seqnn
    from baskerville import dna
except ImportError as e:
    print(f"[FATAL] Baskerville module not found. Ensure conda env is active: {e}")
    sys.exit(1)

# ------------------------------------------------------------------------------
# User Configuration
# ------------------------------------------------------------------------------
WT_FASTA     = "WT_Sniffles2.DEL.2E9M11_524k.fa"   
MUT_FASTA    = "MUT_Sniffles2.DEL.2E9M11_524k.fa"  
PARAMS_FILE  = "params_pred.json"
MODEL_FILE   = "model0_best.h5"
GFF_FILE     = "GCF_000003025.6_Sscrofa11.1_genomic.gff" 

TARGET_GENE  = "CNOT2"
CHROM        = "NC_010447.5"
SV_POS       = 34443468    # Absolute physical center coordinate of SV
SV_TYPE      = "Deletion"
SV_LEN_KB    = 0.294       
HALF_SV_KB   = SV_LEN_KB / 2.0
VIEW_WINDOW_KB = 250       # Window size for plotting (±250 kb)

PLOT_CAGE    = True        

# ------------------------------------------------------------------------------
# Visualization Utilities
# ------------------------------------------------------------------------------
def read_fasta(fasta_path):
    with open(fasta_path, 'r') as f:
        lines = f.read().strip().split('\n')
        return "".join(lines[1:]).upper()

def apply_nature_theme(ax):
    """Applies a minimalist, high-contrast theme."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    ax.tick_params(axis='both', which='major', labelsize=12, width=1.2, length=6)
    ax.yaxis.grid(True, linestyle='--', color='gray', alpha=0.3)
    ax.xaxis.grid(False)

# ------------------------------------------------------------------------------
# GFF Parsing & Gene Track Rendering
# ------------------------------------------------------------------------------
def extract_genes_from_gff(chrom, start, end):
    """Dynamically parses GFF to extract coding elements within the plotting window."""
    print(f"[*] Scanning GFF database for elements in {chrom}:{start}-{end}...")
    genes = {}
    
    with open(GFF_FILE, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            if parts[0] != chrom: continue
            
            f_start, f_end = int(parts[3]), int(parts[4])
            if f_end < start or f_start > end: continue
            
            f_type, strand, info = parts[2], parts[6], parts[8]
            
            if f_type == "gene":
                gene_name = "Unknown"
                for item in info.split(';'):
                    if item.startswith("Name="):
                        gene_name = item.split('=')[1]
                if gene_name not in genes:
                    genes[gene_name] = {"name": gene_name, "strand": strand, "start": f_start, "end": f_end, "exons": []}
            elif f_type == "exon":
                for item in info.split(';'):
                    if item.startswith("gene="):
                        g_name = item.split('=')[1]
                        if g_name in genes:
                            genes[g_name]["exons"].append((f_start, f_end))
                            
    for g in genes.values():
        if not g["exons"]: g["exons"].append((g["start"], g["end"]))
    return list(genes.values())

def plot_gene_track(ax, sv_center_pos, window_kb):
    """Renders the annotated genomic architecture aligned to predictive tracks."""
    abs_start = sv_center_pos - int(window_kb * 1000)
    abs_end = sv_center_pos + int(window_kb * 1000)
    genes = extract_genes_from_gff(CHROM, abs_start, abs_end)
    
    to_rel_kb = lambda pos: (pos - sv_center_pos) / 1000.0
    y_step = -0.6
    
    for i, g in enumerate(genes):
        y_pos = (i % 3) * y_step 
        start_kb = max(-window_kb, to_rel_kb(g["start"]))
        end_kb = min(window_kb, to_rel_kb(g["end"]))
        direction = "<" if g["strand"] == "-" else ">"
        
        # Backbone line & Directional arrows
        ax.plot([start_kb, end_kb], [y_pos, y_pos], color="#1B4F72", lw=1.5)
        arrow_pts = np.linspace(start_kb, end_kb, num=max(2, int((end_kb - start_kb) / 20)))
        for pt in arrow_pts:
            if -window_kb < pt < window_kb:
                ax.text(pt, y_pos, direction, color="#1B4F72", fontsize=10, ha='center', va='center', fontweight='bold')

        # Exon blocks
        for (e_start, e_end) in g["exons"]:
            ex_st_kb = to_rel_kb(e_start)
            ex_width = max(0.2, to_rel_kb(e_end) - ex_st_kb)
            if ex_st_kb + ex_width > -window_kb and ex_st_kb < window_kb:
                ax.add_patch(plt.Rectangle((ex_st_kb, y_pos - 0.1), ex_width, 0.2, lw=0, facecolor='#1B4F72'))
            
        # Gene Label (Highlight Target)
        lbl_color = "#D32F2F" if g['name'] == TARGET_GENE else "black"
        ax.text((start_kb + end_kb)/2, y_pos + 0.15, g['name'], fontsize=11, fontstyle='italic', fontweight='bold', color=lbl_color, ha='center')

    ax.set_ylim(-1.8, 0.5)
    ax.set_xlim(-window_kb, window_kb)
    ax.set_yticks([]); ax.xaxis.grid(False); ax.yaxis.grid(False)
    for spine in ['top', 'right', 'left', 'bottom']: ax.spines[spine].set_visible(False)

# ------------------------------------------------------------------------------
# Main Inference Logic
# ------------------------------------------------------------------------------
def main():
    print("[1] Initializing Deep Learning Architecture...")
    with open(PARAMS_FILE) as f: params = json.load(f)
    model = seqnn.SeqNN(params['model']); model.restore(MODEL_FILE)

    print("[2] Encoding sequences...")
    wt_1hot = np.expand_dims(dna.dna_1hot(read_fasta(WT_FASTA)), axis=0)
    mut_1hot = np.expand_dims(dna.dna_1hot(read_fasta(MUT_FASTA)), axis=0)

    print("[3] Executing In Silico Mutagenesis (Forward Pass)...")
    wt_pred, mut_pred = model.predict(wt_1hot)[0], model.predict(mut_1hot)[0]

    # Target indices (Liver tracks for demonstration)
    rna_p, rna_m = 6491, 6490
    atac_p, atac_m = 1724, None 
    cage_p, cage_m = 23, 22

    print(f"[4] Rendering Epigenomic Projection Landscape (±{VIEW_WINDOW_KB} kb)...")
    center = 16384 // 2
    window = int((VIEW_WINDOW_KB * 1000) / 32)
    start_bin, end_bin = center - window, center + window
    x_axis_kb = np.arange(-window * 32, window * 32, 32) / 1000.0

    n_plots = 2 + (1 if PLOT_CAGE else 0) + 1 
    fig, axes = plt.subplots(n_plots, 1, figsize=(14, 3.5 * n_plots), sharex=True, gridspec_kw={'height_ratios': [3, 3] + ([3] if PLOT_CAGE else []) + [1]})

    plot_idx = 0
    colors = {"wt": "gray", "rna": "#C0392B", "atac": "#27AE60", "cage": "#2980B9", "shade": "#FADBD8"}

    def add_track(ax, p_idx, m_idx, color_key, title, is_rna=False):
        apply_nature_theme(ax)
        ax.plot(x_axis_kb, wt_pred[start_bin:end_bin, p_idx], color=colors["wt"], lw=1.2, label="Reference (WT) (+)")
        ax.plot(x_axis_kb, mut_pred[start_bin:end_bin, p_idx], color=colors[color_key], lw=1.2, alpha=0.9, label=f"Mutation ({SV_LEN_KB}kb {SV_TYPE}) (+)")
        if m_idx is not None:
            ax.plot(x_axis_kb, wt_pred[start_bin:end_bin, m_idx], color=colors["wt"], lw=1.2, ls='--', alpha=0.7, label="Reference (WT) (-)")
            ax.plot(x_axis_kb, mut_pred[start_bin:end_bin, m_idx], color=colors[color_key], lw=1.2, ls='--', alpha=0.7, label="Mutation (-)")
        ax.set_title(title, fontsize=15, fontweight='bold', pad=15)
        ax.set_ylabel("Predicted Signal", fontsize=13, fontweight='bold')
        if is_rna: ax.legend(loc="upper right", frameon=False, fontsize=11, ncol=2)
        ax.axvspan(-HALF_SV_KB, HALF_SV_KB, color=colors["shade"], alpha=0.6, zorder=0)

    add_track(axes[plot_idx], rna_p, rna_m, "rna", "Predicted RNA-seq Output in Liver", is_rna=True); plot_idx += 1
    add_track(axes[plot_idx], atac_p, atac_m, "atac", "Predicted Chromatin Accessibility", is_rna=False); plot_idx += 1
    if PLOT_CAGE:
        add_track(axes[plot_idx], cage_p, cage_m, "cage", "Predicted Promoter Activity (CAGE)", is_rna=False); plot_idx += 1

    ax_g = axes[plot_idx]
    plot_gene_track(ax_g, SV_POS, VIEW_WINDOW_KB)
    ax_g.set_xlabel(f"Genomic Distance from target site (kb)", fontsize=15, fontweight='bold', labelpad=15)
    
    for a in axes: a.axvline(0, color='red', ls='--', lw=1.5, alpha=0.7, zorder=1)

    plt.subplots_adjust(hspace=0.3) 
    out_pdf = f"InSilico_Mutagenesis_Landscape_{TARGET_GENE}.pdf"
    plt.savefig(out_pdf, bbox_inches='tight')
    print(f"[+] Success! High-resolution vector plot saved to: {out_pdf}")

if __name__ == '__main__':
    main()
