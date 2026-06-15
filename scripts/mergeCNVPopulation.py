#! /usr/bin/env python
"""Cross-sample CNV merge -> population CN / binary matrices (standalone CLI).

Split out of the Snakemake workflow so each species can be merged with its own
overlap parameters, typically after per-sample refiltering with filterCNV.py.

Reuses the same merge core as the old in-workflow populationMatrix.py
(utils.mergeCNVFromTools + calculateOverlapProp4Region). Defaults match the
original config.yaml (min-overlap 0.75, max-overlap 0.95, length-ratio 0.5,
freq-threshold 0.8).

Note: mergeCNVFromTools is a greedy front-to-back merge, so region boundaries
depend on the order CNVs are fed in. By default samples are processed in sorted
order; pass --samples in the original config order to reproduce a prior matrix
byte-for-byte.

Input: per-sample BEDs (default res/filtered/). The first 4 columns
(chromosome, start, end, cn) are used; the header line is skipped.

Output:
    <out-dir>/cn_matrix.tsv       per-region copy number per sample
    <out-dir>/binary_matrix.tsv   per-region presence/absence per sample
Both carry columns: chrom start end cnv_type pop_freq ref_bias_flag <samples...>
"""
import argparse
import sys
from pathlib import Path

# Resolve FissionCNV/scripts (this file's directory) for utils import.
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from utils import mergeCNVFromTools, calculateOverlapProp4Region, NORMAL_CN

# Defaults mirror the original config.yaml merge block.
DEF_MIN_OVERLAP = 0.75
DEF_MAX_OVERLAP = 0.95
DEF_LENGTH_RATIO = 0.5
DEF_FREQ_THRESHOLD = 0.8


def load_per_sample_cnvs(in_dir, sample_names):
    """Return (sample -> [[chrom,start,end,cn],...], flat list with sample col)."""
    sample_cnvs = {}
    all_cnvs = []
    for s in sample_names:
        p = in_dir / f"{s}.bed"
        rows = []
        if p.exists():
            with open(p) as f:
                for line in f:
                    if line.startswith("chrom"):
                        continue
                    x = line.rstrip("\n").split("\t")
                    if len(x) < 4:
                        continue
                    chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
                    rows.append([chrom, start, end, cn])
                    all_cnvs.append([chrom, start, end, cn, s])
        else:
            print(f"[warn] missing {p}", flush=True)
        sample_cnvs[s] = rows
    return sample_cnvs, all_cnvs


def build_matrices(sample_cnvs, all_cnvs, sample_names, out_dir,
                   min_overlap, max_overlap, length_ratio, freq_threshold):
    cnv_for_merge = [[c[0], c[1], c[2], c[3], c[4]] for c in all_cnvs]
    unified = mergeCNVFromTools(
        cnv_for_merge,
        min_threshold=min_overlap,
        max_threshold=max_overlap,
        length_ratio_limit=length_ratio,
    )
    # unified: [chrom, start, end, cn, samples_str, sample_count, accumScore]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cn = out_dir / "cn_matrix.tsv"
    out_bin = out_dir / "binary_matrix.tsv"

    n = len(sample_names)
    header = ["chrom", "start", "end", "cnv_type", "pop_freq", "ref_bias_flag"] + sample_names
    with open(out_cn, "w") as fc, open(out_bin, "w") as fb:
        fc.write("\t".join(header) + "\n")
        fb.write("\t".join(header) + "\n")
        for region in unified:
            chrom, start, end, cn = region[0], region[1], region[2], region[3]
            cnv_type = "DEL" if cn < NORMAL_CN else "DUP"
            region_coord = [chrom, start, end]
            cn_row, bin_row = [], []
            for s in sample_names:
                found_cn, found = NORMAL_CN, 0
                for c in sample_cnvs[s]:
                    if c[0] != chrom:
                        continue
                    overlap, p1, p2 = calculateOverlapProp4Region(region_coord, [c[0], c[1], c[2]])
                    if min(p1, p2) > min_overlap or max(p1, p2) > max_overlap:
                        if (c[3] < NORMAL_CN and cn < NORMAL_CN) or \
                           (c[3] > NORMAL_CN and cn > NORMAL_CN):
                            found_cn, found = c[3], 1
                            break
                cn_row.append(str(found_cn))
                bin_row.append(str(found))
            pop_freq = round(sum(int(x) for x in bin_row) / n, 4) if n else 0.0
            ref_bias = "YES" if pop_freq > freq_threshold else "NO"
            prefix = [chrom, str(start), str(end), cnv_type, str(pop_freq), ref_bias]
            fc.write("\t".join(prefix + cn_row) + "\n")
            fb.write("\t".join(prefix + bin_row) + "\n")
    return out_cn, out_bin, len(unified)


def main():
    ap = argparse.ArgumentParser(
        description="Cross-sample CNV merge into population CN/binary matrices.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--in-dir", default="res/filtered",
                    help="Directory of per-sample BEDs to merge")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for the two matrices")
    ap.add_argument("--samples", nargs="*",
                    help="Sample names (without .bed); default: all *.bed in in-dir")
    ap.add_argument("--min-overlap", type=float, default=DEF_MIN_OVERLAP,
                    help="Minimum reciprocal overlap to merge two regions")
    ap.add_argument("--max-overlap", type=float, default=DEF_MAX_OVERLAP,
                    help="One-sided overlap above which regions merge regardless")
    ap.add_argument("--length-ratio", type=float, default=DEF_LENGTH_RATIO,
                    help="Min(len)/max(len) below which a pair is not merged")
    ap.add_argument("--freq-threshold", type=float, default=DEF_FREQ_THRESHOLD,
                    help="pop_freq above this flags ref_bias_flag=YES")
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    if not in_dir.is_dir():
        ap.error(f"--in-dir not found: {in_dir}")

    if args.samples:
        sample_names = args.samples
    else:
        sample_names = sorted(p.stem for p in in_dir.glob("*.bed"))
    if not sample_names:
        ap.error(f"no input BEDs found in {in_dir}")

    print(f"[info] {len(sample_names)} samples from {in_dir}", flush=True)
    print(f"[info] min_overlap={args.min_overlap} max_overlap={args.max_overlap} "
          f"length_ratio={args.length_ratio} freq_threshold={args.freq_threshold}",
          flush=True)

    sample_cnvs, all_cnvs = load_per_sample_cnvs(in_dir, sample_names)
    print(f"[info] cross-sample merge of {len(all_cnvs)} CNV records", flush=True)
    out_cn, out_bin, n_regions = build_matrices(
        sample_cnvs, all_cnvs, sample_names, out_dir,
        args.min_overlap, args.max_overlap, args.length_ratio, args.freq_threshold)
    print(f"[done] {n_regions} unified regions", flush=True)
    print(f"[done] {out_cn}", flush=True)
    print(f"[done] {out_bin}", flush=True)


if __name__ == "__main__":
    main()
