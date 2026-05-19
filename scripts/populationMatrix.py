#! /usr/bin/env python
# Generate population CNV matrix from per-sample final BED files.
# Output: CN matrix and binary (presence/absence) matrix.

import sys
import os
sys.path.insert(0, snakemake.params['absPath'] + '/scripts')
from utils import mergeCNVFromTools, calculateOverlapProp4Region, NORMAL_CN

sample_names = snakemake.params['sample_names']
min_overlap = snakemake.params['min_overlap']
max_overlap = snakemake.params['max_overlap']
length_ratio = snakemake.params['length_ratio']
freq_threshold = snakemake.params['freq_threshold']
output_cn = snakemake.output['cn_matrix']
output_binary = snakemake.output['binary_matrix']

# Step 1: Read all per-sample CNVs
all_cnvs = []  # list of [chrom, start, end, cn, sample_name]
sample_cnvs = {}  # sample -> list of [chrom, start, end, cn]

for sample in sample_names:
    bed_file = f'res/final/{sample}.bed'
    sample_cnvs[sample] = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith('chromosome'):
                continue
            x = line.strip().split('\t')
            if len(x) < 4:
                continue
            chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
            sample_cnvs[sample].append([chrom, start, end, cn])
            all_cnvs.append([chrom, start, end, cn, sample])

# Step 2: Cross-sample merge to define unified CNV regions
# Prepare for merging: add sample as tool name (reuse mergeCNVFromTools)
cnv_for_merge = [[c[0], c[1], c[2], c[3], c[4]] for c in all_cnvs]
unified_regions = mergeCNVFromTools(
    cnv_for_merge,
    min_threshold=min_overlap,
    max_threshold=max_overlap,
    length_ratio_limit=length_ratio
)
# unified_regions: [chrom, start, end, cn, samples_str, sample_count, accumScore]

# Step 3: Build matrices
# For each unified region, check each sample for overlap
os.makedirs(os.path.dirname(output_cn), exist_ok=True)

# Build CN matrix and binary matrix
with open(output_cn, 'w') as f_cn, open(output_binary, 'w') as f_bin:
    # Header
    header = ['chrom', 'start', 'end', 'cnv_type', 'pop_freq', 'ref_bias_flag'] + sample_names
    print('\t'.join(header), file=f_cn)
    print('\t'.join(header), file=f_bin)

    for region in unified_regions:
        chrom, start, end, cn = region[0], region[1], region[2], region[3]
        cnv_type = 'DEL' if cn < NORMAL_CN else 'DUP'
        region_coord = [chrom, start, end]

        cn_row = []
        bin_row = []

        for sample in sample_names:
            found_cn = NORMAL_CN  # default: normal
            found = 0
            for s_cnv in sample_cnvs[sample]:
                # Check overlap between sample CNV and unified region
                s_coord = [s_cnv[0], s_cnv[1], s_cnv[2]]
                overlap, prop1, prop2 = calculateOverlapProp4Region(region_coord, s_coord)
                if min(prop1, prop2) > min_overlap or max(prop1, prop2) > max_overlap:
                    # Check same CNV type
                    if (s_cnv[3] < NORMAL_CN and cn < NORMAL_CN) or \
                       (s_cnv[3] > NORMAL_CN and cn > NORMAL_CN):
                        found_cn = s_cnv[3]
                        found = 1
                        break
            cn_row.append(str(found_cn))
            bin_row.append(str(found))

        # Calculate population frequency
        pop_freq = round(sum(int(x) for x in bin_row) / len(sample_names), 4)
        ref_bias = 'YES' if pop_freq > freq_threshold else 'NO'

        row_prefix = [chrom, str(start), str(end), cnv_type, str(pop_freq), ref_bias]
        print('\t'.join(row_prefix + cn_row), file=f_cn)
        print('\t'.join(row_prefix + bin_row), file=f_bin)
