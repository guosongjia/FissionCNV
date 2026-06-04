#!/usr/bin/env python
"""Unified breakpoint annotation: appends breakpoint context columns to each CNV.

Precise-tool CNVs (smoove/delly): CI extraction + microhomology + repeat + BAM evidence.
Depth-only CNVs (cnvkit/cnvpytor/mops): repeat + BAM evidence only (mh unreliable at bin-level).
"""

import sys
import subprocess
import tempfile
import os
import shutil
import pandas as pd
from pyfaidx import Fasta
import pysam

sys.path.insert(0, snakemake.params['absPath'] + '/scripts')

DEPTH_TOOLS = {'cnvkit', 'cnvpytor', 'mops'}
PRECISE_TOOLS = {'smoove', 'delly'}
MIN_WINDOW = 10
MAX_MH = 25
REPEAT_PRIORITY = ['KMDs', 'LTR', 'tRNA', 'rRNA', 'centromeric', 'subtelomeric']


def parse_ci_from_vcf(vcf_path, tool):
    ci_map = {}
    fmt = '%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/CIPOS\\t%INFO/CIEND\\n'
    cmd = ['bcftools', 'query', '-f', fmt, vcf_path]
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ci_map
    for line in out.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 6:
            continue
        chrom, pos, end, svtype, cipos_str, ciend_str = parts[:6]
        if svtype not in ('DEL', 'DUP'):
            continue
        cipos = _parse_ci_field(cipos_str)
        ciend = _parse_ci_field(ciend_str)
        ci_map[(chrom, int(pos) - 1, int(end), svtype)] = (cipos, ciend)
    return ci_map


def _parse_ci_field(s):
    if not s or s == '.':
        return (0, 0)
    parts = s.split(',')
    if len(parts) == 2:
        return (abs(int(parts[0])), abs(int(parts[1])))
    return (0, 0)

def find_ci_for_cnv(chrom, start, end, svtype, ci_maps):
    best_cipos = (0, 0)
    best_ciend = (0, 0)
    best_dist = float('inf')
    for ci_map in ci_maps:
        for (c, s, e, st), (cipos, ciend) in ci_map.items():
            if c != chrom or st != svtype:
                continue
            dist = abs(s - start) + abs(e - end)
            if dist < best_dist:
                best_dist = dist
                best_cipos = cipos
                best_ciend = ciend
    return best_cipos, best_ciend


def compute_microhomology(ref, chrom, start, end, svtype):
    try:
        seq = ref[chrom]
    except KeyError:
        return 0, 0, ''
    chrom_len = len(seq)
    if svtype == 'DUP':
        left_pos, right_pos = end, start
    else:
        left_pos, right_pos = start, end

    strict_mh = 0
    fuzzy_mh = 0
    mh_seq_chars = []
    mismatches = 0

    for k in range(MAX_MH):
        lp = left_pos + k
        rp = right_pos + k
        if lp >= chrom_len or rp >= chrom_len:
            break
        base_l = str(seq[lp]).upper()
        base_r = str(seq[rp]).upper()
        if base_l == base_r:
            if mismatches == 0:
                strict_mh = k + 1
            fuzzy_mh = k + 1
            mh_seq_chars.append(base_l)
        else:
            mismatches += 1
            if mismatches > 1:
                break
            fuzzy_mh = k + 1
            mh_seq_chars.append(base_l.lower())

    mh_seq = ''.join(mh_seq_chars[:strict_mh]) if strict_mh > 0 else ''
    return strict_mh, fuzzy_mh, mh_seq


def intersect_repeats(bp_bed_path, annotation_files):
    hits = {}
    for cat, path in annotation_files.items():
        if not path or not os.path.isfile(path):
            continue
        cmd = ['bedtools', 'intersect', '-a', bp_bed_path, '-b', path, '-wa', '-wb']
        try:
            out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True)
        except subprocess.CalledProcessError:
            continue
        for line in out.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            interval_key = (parts[0], int(parts[1]), int(parts[2]), parts[3])
            if interval_key not in hits:
                hits[interval_key] = []
            hits[interval_key].append(cat)
    return hits


def get_priority_hit(hit_list):
    for cat in REPEAT_PRIORITY:
        if cat in hit_list:
            return cat
    return 'none'


def check_nested(cnvs, idx):
    target = cnvs[idx]
    for i, other in enumerate(cnvs):
        if i == idx:
            continue
        if other[0] != target[0]:
            continue
        if other[1] <= target[1] and other[2] >= target[2] and (other[1] < target[1] or other[2] > target[2]):
            return True
    return False


def scan_bam_evidence(bam_path, chrom, pos, window):
    split_reads = 0
    discordant_pairs = 0
    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            for read in bam.fetch(chrom, max(0, pos - window), pos + window):
                if read.is_unmapped or read.is_secondary or read.is_duplicate:
                    continue
                if read.has_tag('SA'):
                    split_reads += 1
                elif read.is_paired and not read.is_proper_pair and not read.mate_is_unmapped:
                    discordant_pairs += 1
    except (ValueError, OSError):
        pass
    return split_reads, discordant_pairs


# === Main ===
sample = snakemake.wildcards.sample
absPath = snakemake.params['absPath']
ref_path = snakemake.params['genome']
bam_path = snakemake.input['bam']
input_bed = snakemake.input['annotated_bed']
output_bed = snakemake.output[0]

annotation_files = {}
ann_map = {
    'trna_gff': 'tRNA', 'rrna_gff': 'rRNA', 'ltr_bed': 'LTR',
    'centromeric_bed': 'centromeric', 'subtelomeric_bed': 'subtelomeric', 'kmds_bed': 'KMDs',
}
for key, cat in ann_map.items():
    val = snakemake.params[key]
    if val:
        annotation_files[cat] = val

longdust_path = snakemake.params['low_complexity']

ref = Fasta(ref_path)
chrom_lengths = {name: len(ref[name]) for name in ref.keys()}

# Load raw VCF CI data for precise tools
ci_maps = []
smoove_vcf = snakemake.input.get('smoove_vcf', '')
delly_bcf = snakemake.input.get('delly_bcf', '')
if smoove_vcf and os.path.isfile(smoove_vcf):
    ci_maps.append(parse_ci_from_vcf(smoove_vcf, 'smoove'))
if delly_bcf and os.path.isfile(delly_bcf):
    ci_maps.append(parse_ci_from_vcf(delly_bcf, 'delly'))

# Read annotated BED — process ALL CNVs
cnvs = []
original_lines = []
with open(input_bed) as f:
    header = f.readline().strip()
    for line in f:
        original_lines.append(line.strip())
        x = line.strip().split('\t')
        chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
        svtype = 'DEL' if cn < 1 else 'DUP'
        tool_str = x[10]
        tool_num = int(x[11])
        has_precise = bool(set(tool_str.split(',')) & PRECISE_TOOLS)
        cnvs.append((chrom, start, end, cn, svtype, tool_str, tool_num, has_precise))

# Build breakpoint windows for repeat intersection
tmpdir = tempfile.mkdtemp()
bp_bed_path = os.path.join(tmpdir, 'bp_windows.bed')
with open(bp_bed_path, 'w') as bp_f:
    for i, (chrom, start, end, cn, svtype, tools, tnum, has_precise) in enumerate(cnvs):
        if has_precise:
            cipos, ciend = find_ci_for_cnv(chrom, start, end, svtype, ci_maps)
        else:
            cipos, ciend = (0, 0), (0, 0)
        w5_left = max(0, start - max(cipos[0], MIN_WINDOW))
        w5_right = min(chrom_lengths.get(chrom, end), start + max(cipos[1], MIN_WINDOW))
        w3_left = max(0, end - max(ciend[0], MIN_WINDOW))
        w3_right = min(chrom_lengths.get(chrom, end), end + max(ciend[1], MIN_WINDOW))
        bp_f.write(f"{chrom}\t{w5_left}\t{w5_right}\t{i}_5p\n")
        bp_f.write(f"{chrom}\t{w3_left}\t{w3_right}\t{i}_3p\n")

repeat_hits = intersect_repeats(bp_bed_path, annotation_files)
simple_hits = {}
if longdust_path and os.path.isfile(longdust_path):
    simple_hits = intersect_repeats(bp_bed_path, {'simple': longdust_path})

# Per-CNV analysis: build breakpoint columns
bp_columns = []
for i, (chrom, start, end, cn, svtype, tools, tnum, has_precise) in enumerate(cnvs):
    if has_precise:
        cipos, ciend = find_ci_for_cnv(chrom, start, end, svtype, ci_maps)
    else:
        cipos, ciend = (0, 0), (0, 0)
    w5_left = max(0, start - max(cipos[0], MIN_WINDOW))
    w5_right = min(chrom_lengths.get(chrom, end), start + max(cipos[1], MIN_WINDOW))
    w3_left = max(0, end - max(ciend[0], MIN_WINDOW))
    w3_right = min(chrom_lengths.get(chrom, end), end + max(ciend[1], MIN_WINDOW))

    at_chrom_boundary = (w5_left == 0) or (w3_right == chrom_lengths.get(chrom, end))
    nested = check_nested(cnvs, i)

    key_5p = (chrom, w5_left, w5_right, f"{i}_5p")
    key_3p = (chrom, w3_left, w3_right, f"{i}_3p")
    hits_5p = repeat_hits.get(key_5p, [])
    hits_3p = repeat_hits.get(key_3p, [])
    repeat_type_5p = get_priority_hit(hits_5p)
    repeat_type_3p = get_priority_hit(hits_3p)
    repeat_types_all_5p = ','.join(sorted(set(hits_5p))) if hits_5p else 'none'
    repeat_types_all_3p = ','.join(sorted(set(hits_3p))) if hits_3p else 'none'
    simple_5p = 'simple' in simple_hits.get(key_5p, [])
    simple_3p = 'simple' in simple_hits.get(key_3p, [])
    in_repeat = bool(hits_5p) or bool(hits_3p)

    # BAM evidence (always)
    sr_5p, disc_5p = scan_bam_evidence(bam_path, chrom, start, MIN_WINDOW)
    sr_3p, disc_3p = scan_bam_evidence(bam_path, chrom, end, MIN_WINDOW)

    # Microhomology (precise tools only)
    if has_precise:
        length = end - start
        dup_too_short = (svtype == 'DUP' and length < 50)
        if dup_too_short:
            strict_mh, fuzzy_mh, mh_seq = 0, 0, ''
        else:
            strict_mh, fuzzy_mh, mh_seq = compute_microhomology(ref, chrom, start, end, svtype)
        dup_assumption = 'tandem' if svtype == 'DUP' else ''
    else:
        strict_mh, fuzzy_mh, mh_seq = '', '', ''
        dup_assumption = ''
        dup_too_short = ''

    bp_columns.append('\t'.join(str(v) for v in [
        cipos[0], cipos[1], ciend[0], ciend[1],
        strict_mh, fuzzy_mh, mh_seq, dup_assumption,
        in_repeat, repeat_type_5p, repeat_type_3p,
        repeat_types_all_5p, repeat_types_all_3p,
        simple_5p, simple_3p,
        sr_5p, sr_3p, disc_5p, disc_3p,
        nested, at_chrom_boundary, dup_too_short,
    ]))

# Write output: original lines + breakpoint columns
BP_HEADER = '\t'.join([
    'ci_5p_left', 'ci_5p_right', 'ci_3p_left', 'ci_3p_right',
    'strict_mh_len', 'fuzzy_mh_len', 'mh_seq', 'dup_assumption',
    'in_repeat', 'repeat_type_5p', 'repeat_type_3p',
    'repeat_types_all_5p', 'repeat_types_all_3p',
    'simple_overlap_5p', 'simple_overlap_3p',
    'split_reads_5p', 'split_reads_3p', 'discordant_pairs_5p', 'discordant_pairs_3p',
    'nested', 'at_chrom_boundary', 'dup_too_short',
])

with open(output_bed, 'w') as f:
    f.write(header + '\t' + BP_HEADER + '\n')
    for orig_line, bp_cols in zip(original_lines, bp_columns):
        f.write(orig_line + '\t' + bp_cols + '\n')

shutil.rmtree(tmpdir, ignore_errors=True)
