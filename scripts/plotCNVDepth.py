#!/usr/bin/env python
"""
Standalone script to generate per-CNV read depth plots for a specified strain.

Usage:
    python scripts/plotCNVDepth.py --bam mapped/strain.bam \
        --bed res/final/strain.bed --fai genome.fa.fai --output plots/

    # Plot a single region:
    python scripts/plotCNVDepth.py --bam mapped/strain.bam \
        --region chr1:1000-5000 --fai genome.fa.fai --output plots/
"""

import argparse
import subprocess
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def get_depth(bam, chrom, start, end):
    """Get per-base depth using samtools depth."""
    cmd = f"samtools depth -r {chrom}:{start}-{end} {bam}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    positions = []
    depths = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        positions.append(int(parts[1]))
        depths.append(int(parts[2]))
    return positions, depths


def plot_cnv_depth(bam, chrom, start, end, cn, output_dir, fai_file):
    """Plot read depth for a single CNV region with flanking."""
    cnv_len = end - start
    flank = max(int(cnv_len * 0.3), 500)
    plot_start = max(0, start - flank)
    plot_end = end + flank

    positions, depths = get_depth(bam, chrom, plot_start, plot_end)

    if not positions:
        return

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.fill_between(positions, depths, alpha=0.4, color='steelblue')
    ax.axvspan(start, end, alpha=0.2, color='red' if cn == 0 else 'blue')
    ax.axvline(start, color='black', linewidth=0.5, linestyle='--')
    ax.axvline(end, color='black', linewidth=0.5, linestyle='--')

    ax.set_xlabel('Position')
    ax.set_ylabel('Read Depth')
    cnv_type = 'DEL' if cn < 1 else 'DUP'
    ax.set_title(f'{chrom}:{start}-{end} (CN={cn}, {cnv_type})')

    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f'{chrom}_{start}_{end}_{cnv_type}.png')
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close()
    print(f"  Saved: {out_file}")


def main():
    parser = argparse.ArgumentParser(description='Plot per-CNV read depth')
    parser.add_argument('--bam', required=True, help='BAM file path')
    parser.add_argument('--bed', help='CNV BED file (plot all CNVs)')
    parser.add_argument('--region', help='Single region (chr:start-end)')
    parser.add_argument('--fai', required=True, help='Genome .fai file')
    parser.add_argument('--output', required=True, help='Output directory')
    args = parser.parse_args()

    if args.region:
        chrom, coords = args.region.split(':')
        start, end = [int(x) for x in coords.split('-')]
        plot_cnv_depth(args.bam, chrom, start, end, 0, args.output, args.fai)
    elif args.bed:
        with open(args.bed) as f:
            for line in f:
                if line.startswith('chromosome'):
                    continue
                x = line.strip().split('\t')
                if len(x) < 4:
                    continue
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
                plot_cnv_depth(args.bam, chrom, start, end, cn, args.output, args.fai)
    else:
        parser.error("Either --bed or --region must be provided")


if __name__ == '__main__':
    main()
