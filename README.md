# FissionCNV

A Snakemake pipeline for CNV calling in haploid fission yeast (*Schizosaccharomyces*) populations. Adapted from [CNVPipe](https://github.com/sunjh22/CNVPipe).

## Workflow

1. **Preprocessing** (optional): fastp → BWA MEM → markdup
2. **CNV calling** with 5 tools in parallel:
   - [CNVKit](https://github.com/etal/cnvkit) — flat-reference mode, ploidy 1
   - [CNVpytor](https://github.com/abyzovlab/CNVpytor) — read-depth, ploidy 1
   - [cn.MOPS](http://bioconductor.org/packages/cn.mops/) — `haplocn.mops()`
   - [Smoove](https://github.com/brentp/smoove) — split-read SV
   - [Delly](https://github.com/dellytools/delly) — split-read + read-pair SV
3. **Merge** calls across tools (priority: smoove → delly → cnvkit → cnvpytor → mops)
4. **Score** with [duphold](https://github.com/brentp/duphold)
5. **Hard filter** (toolNum ≥ 2, duphold > 0, depth-only ≥ 600 bp)
6. **Annotate** with low-complexity (GS) and low-mappability (MS) overlap, plus population frequency
7. **Output**: per-sample BED, population CN matrix, per-sample genome overview plot

## Usage

```bash
# 1. Edit config.yaml: set genome, samples.tsv path, annotation BEDs, outdir
# 2. Edit samples.tsv: list sample names and fastq paths
# 3. Run
snakemake --profile profiles/default
```

The profile sets `--use-conda`, `--conda-prefix conda_envs`, `--latency-wait 60`, `-j 10`.

## Key config options

```yaml
data:
  samples: samples.tsv
  genome: /path/to/genome.fasta
  outdir: .                       # all output goes here
  low-complexity: /path/to/longdust.bed
  low-mappable: /path/to/genmap.bed

params:
  bam-input: False                # set True to skip preprocessing
  skip-fastp: False               # set True if input reads are pre-cleaned
  binSize: 200
  ploidy: 1
```

## Differences from CNVPipe

- Haploid: CN=0 (DEL), CN=1 (normal), CN≥2 (DUP). Smoove/Delly diploid CN is halved.
- No single-cell, no BAF filter, no ClassifyCNV, no GATK/freebayes, no autobin.
- Each species runs independently with its own reference and config.

## License

Inherits from upstream CNVPipe.
