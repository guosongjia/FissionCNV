#! /usr/bin/env python
"""Per-sample CNV refilter on duphold fold-change metrics (standalone CLI).

Decoupled from the Snakemake workflow so different species / experiments can be
refiltered reproducibly with different thresholds before cross-sample merging
(see mergeCNVPopulation.py).

Input BEDs are the per-sample, already hard-filtered calls in `res/filtered/`
(toolNum >= min-tools, dupholdScore > 0, depth-only >= min-length). The annotated
(`res/annotated/`) and final (`res/final/`) BEDs share the same leading columns
and also work as input.

Relevant columns (0-based):
    4  cnv          "duplication" | "deletion"
    7  dhfc         duphold fold-change vs whole chromosome
    8  dhbfc        duphold fold-change vs GC-matched bins
    9  dhffc        duphold fold-change vs flanking regions

Keep rule (per metric, only applied when a threshold is given):
    DUP keeps a call when its metric is ABOVE the DUP threshold.
    DEL keeps a call when its metric is BELOW the DEL threshold.
With --mode all (default) every supplied threshold must pass; with --mode any a
single passing threshold is enough. If no threshold is supplied for a type, all
calls of that type pass through unchanged.

Examples
--------
# Tier-1 duphold-strict (DUP all three > 1.3, DEL all three < 0.7):
python scripts/filterCNV.py --in-dir res/filtered --out-dir res/refilt_tier1 \\
    --preset tier1

# Only constrain the GC-corrected metric, DUP only:
python scripts/filterCNV.py --in-dir res/filtered --out-dir res/refilt_dhbfc \\
    --types dup --dup-dhbfc 1.5
"""
import argparse
import json
import sys
from pathlib import Path

# cnv type and the three duphold fold-change columns (0-based)
COL_CNV = 4
COL_DHFC = 7
COL_DHBFC = 8
COL_DHFFC = 9

PRESETS = {
    # name -> (dup thresholds, del thresholds) for (dhfc, dhbfc, dhffc)
    "tier1": dict(dup_dhfc=1.3, dup_dhbfc=1.3, dup_dhffc=1.3,
                  del_dhfc=0.7, del_dhbfc=0.7, del_dhffc=0.7),
}


def build_thresholds(args):
    """Resolve preset first, then let explicit flags override."""
    thr = dict(dup_dhfc=None, dup_dhbfc=None, dup_dhffc=None,
               del_dhfc=None, del_dhbfc=None, del_dhffc=None)
    if args.preset:
        thr.update(PRESETS[args.preset])
    for k in thr:
        v = getattr(args, k)
        if v is not None:
            thr[k] = v
    return thr


def passes(cnv_type, dhfc, dhbfc, dhffc, thr, mode):
    """Return True if a call should be kept under the given thresholds."""
    checks = []
    if cnv_type == "duplication":
        if thr["dup_dhfc"] is not None:
            checks.append(dhfc > thr["dup_dhfc"])
        if thr["dup_dhbfc"] is not None:
            checks.append(dhbfc > thr["dup_dhbfc"])
        if thr["dup_dhffc"] is not None:
            checks.append(dhffc > thr["dup_dhffc"])
    elif cnv_type == "deletion":
        if thr["del_dhfc"] is not None:
            checks.append(dhfc < thr["del_dhfc"])
        if thr["del_dhbfc"] is not None:
            checks.append(dhbfc < thr["del_dhbfc"])
        if thr["del_dhffc"] is not None:
            checks.append(dhffc < thr["del_dhffc"])
    else:
        return False  # unknown type
    if not checks:
        return True  # no constraint for this type -> pass through
    return all(checks) if mode == "all" else any(checks)


def filter_one_bed(in_path, out_path, thr, mode, keep_types):
    """Filter one BED. Returns (n_in, n_kept_dup, n_kept_del, n_dropped_bad)."""
    n_in = n_dup = n_del = n_bad = 0
    with open(in_path) as fi, open(out_path, "w") as fo:
        header = fi.readline()
        fo.write(header)
        for line in fi:
            if not line.strip():
                continue
            n_in += 1
            x = line.rstrip("\n").split("\t")
            cnv_type = x[COL_CNV]
            if cnv_type == "duplication" and "dup" not in keep_types:
                continue
            if cnv_type == "deletion" and "del" not in keep_types:
                continue
            try:
                dhfc = float(x[COL_DHFC])
                dhbfc = float(x[COL_DHBFC])
                dhffc = float(x[COL_DHFFC])
            except (ValueError, IndexError):
                n_bad += 1
                continue
            if passes(cnv_type, dhfc, dhbfc, dhffc, thr, mode):
                fo.write(line)
                if cnv_type == "duplication":
                    n_dup += 1
                else:
                    n_del += 1
    return n_in, n_dup, n_del, n_bad


def main():
    ap = argparse.ArgumentParser(
        description="Refilter per-sample CNV BEDs on duphold fold-change metrics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--in-dir", default="res/filtered",
                    help="Directory of per-sample input BEDs")
    ap.add_argument("--out-dir", required=True,
                    help="Directory for filtered per-sample BEDs")
    ap.add_argument("--samples", nargs="*",
                    help="Sample names (without .bed); default: all *.bed in in-dir")
    ap.add_argument("--types", nargs="+", choices=["dup", "del"],
                    default=["dup", "del"],
                    help="CNV types to keep")
    ap.add_argument("--mode", choices=["all", "any"], default="all",
                    help="Combine multiple thresholds with AND (all) or OR (any)")
    ap.add_argument("--preset", choices=sorted(PRESETS),
                    help="Threshold preset; explicit flags override it")
    # DUP: keep when metric > threshold
    ap.add_argument("--dup-dhfc", type=float, dest="dup_dhfc",
                    help="DUP keep iff dhfc > this")
    ap.add_argument("--dup-dhbfc", type=float, dest="dup_dhbfc",
                    help="DUP keep iff dhbfc > this")
    ap.add_argument("--dup-dhffc", type=float, dest="dup_dhffc",
                    help="DUP keep iff dhffc > this")
    # DEL: keep when metric < threshold
    ap.add_argument("--del-dhfc", type=float, dest="del_dhfc",
                    help="DEL keep iff dhfc < this")
    ap.add_argument("--del-dhbfc", type=float, dest="del_dhbfc",
                    help="DEL keep iff dhbfc < this")
    ap.add_argument("--del-dhffc", type=float, dest="del_dhffc",
                    help="DEL keep iff dhffc < this")
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    if not in_dir.is_dir():
        ap.error(f"--in-dir not found: {in_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    thr = build_thresholds(args)
    keep_types = set(args.types)

    if args.samples:
        sample_names = args.samples
    else:
        sample_names = sorted(p.stem for p in in_dir.glob("*.bed"))
    if not sample_names:
        ap.error(f"no input BEDs found in {in_dir}")

    print(f"[info] thresholds: {thr}", flush=True)
    print(f"[info] mode={args.mode} keep_types={sorted(keep_types)} "
          f"samples={len(sample_names)}", flush=True)

    log_rows = []
    tot_in = tot_dup = tot_del = tot_bad = 0
    for i, s in enumerate(sample_names, 1):
        ip = in_dir / f"{s}.bed"
        if not ip.exists():
            print(f"[warn] missing {ip}", flush=True)
            continue
        op = out_dir / f"{s}.bed"
        n_in, n_dup, n_del, n_bad = filter_one_bed(ip, op, thr, args.mode, keep_types)
        tot_in += n_in; tot_dup += n_dup; tot_del += n_del; tot_bad += n_bad
        log_rows.append((s, n_in, n_dup, n_del, n_dup + n_del, n_bad))
        if i % 100 == 0 or i == len(sample_names):
            print(f"[info] filtered {i}/{len(sample_names)}", flush=True)

    # Reproducibility log: parameters + per-sample counts
    log_path = out_dir / "filter_log.tsv"
    with open(log_path, "w") as f:
        f.write("# filterCNV.py parameters: "
                + json.dumps({"in_dir": str(in_dir), "mode": args.mode,
                              "preset": args.preset, "types": sorted(keep_types),
                              "thresholds": thr}) + "\n")
        f.write("sample\tn_in\tn_kept_dup\tn_kept_del\tn_kept_total\tn_dropped_unparsed\n")
        for r in log_rows:
            f.write("\t".join(str(v) for v in r) + "\n")
        f.write("\t".join(["TOTAL", str(tot_in), str(tot_dup), str(tot_del),
                           str(tot_dup + tot_del), str(tot_bad)]) + "\n")

    print(f"[done] in={tot_in} kept_dup={tot_dup} kept_del={tot_del} "
          f"kept_total={tot_dup + tot_del} dropped_unparsed={tot_bad}", flush=True)
    print(f"[done] per-sample log -> {log_path}", flush=True)


if __name__ == "__main__":
    main()
