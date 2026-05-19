#! /usr/bin/env python
# Generate per-sample whole-genome CNV overview plot.
# X-axis: genomic position (chromosomes as panels)
# Y-axis: copy number
# Colors: CN=0 red (DEL), CN=1 grey (normal), CN>=2 blue (DUP)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

bed_file = snakemake.input[0]
fai_file = snakemake.params['fai']
output_file = snakemake.output[0]
sample = snakemake.wildcards.sample

# Read chromosome lengths from .fai
chrom_lengths = {}
chrom_order = []
with open(fai_file) as f:
    for line in f:
        x = line.strip().split('\t')
        chrom_lengths[x[0]] = int(x[1])
        chrom_order.append(x[0])

# Read CNVs
cnvs = []
with open(bed_file) as f:
    for line in f:
        if line.startswith('chromosome'):
            continue
        x = line.strip().split('\t')
        if len(x) < 4:
            continue
        cnvs.append((x[0], int(x[1]), int(x[2]), int(x[3])))

# Plot
fig, axes = plt.subplots(len(chrom_order), 1, figsize=(12, 2 * len(chrom_order)),
                         squeeze=False, sharex=False)

for idx, chrom in enumerate(chrom_order):
    ax = axes[idx, 0]
    chrom_len = chrom_lengths[chrom]

    # Background: normal CN=1
    ax.axhline(y=1, color='lightgrey', linewidth=0.5, linestyle='--')
    ax.set_xlim(0, chrom_len)
    ax.set_ylim(-0.5, 5.5)
    ax.set_ylabel('CN')
    ax.set_title(f'{chrom} ({chrom_len:,} bp)', fontsize=9, loc='left')

    # Plot CNVs on this chromosome
    for c, s, e, cn in cnvs:
        if c != chrom:
            continue
        color = 'red' if cn == 0 else 'blue'
        ax.axvspan(s, e, alpha=0.6, color=color, linewidth=0)
        ax.plot([s, e], [cn, cn], color=color, linewidth=2)

    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)

plt.suptitle(f'{sample} - CNV Overview', fontsize=12)
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()
