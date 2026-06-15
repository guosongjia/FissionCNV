# FissionCNV

A Snakemake pipeline for CNV calling in haploid fission yeast (*Schizosaccharomyces*) populations. Adapted from [CNVPipe](https://github.com/sunjh22/CNVPipe).

![pipeline overview](pipeline.png)

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
6. **Annotate** with low-complexity (GS) and low-mappability (MS) overlap
7. **Output**: per-sample BED (`res/final/{sample}.bed`), per-sample genome overview plot

Cross-sample population merge is **not** part of the workflow (see below): different
species often need different per-sample refilter and overlap parameters before merging.

### Pipeline diagram

```mermaid
flowchart TD
    %% ====== Inputs ======
    subgraph IN["Input"]
        FQ["paired-end FASTQ<br/>(samples.tsv)"]
        BAM0["existing BAM<br/>(bam-input: True)"]
        REF["reference FASTA<br/>config.data.genome"]
        ANNO["annotation BEDs<br/>low-complexity / low-mappable"]
    end

    %% ====== Preprocessing ======
    subgraph PRE["Preprocessing (skipped when bam-input=True)"]
        FAIDX["samtools faidx<br/>(genome.fai)"]
        BWAIDX["bwa index<br/>(.amb/.ann/.bwt/.pac/.sa)"]
        FASTP["fastp<br/>(skip-fastp: True bypasses)"]
        MULTIQC["multiqc"]
        BWAMEM["bwa mem + samtools fixmate<br/>+ samtools sort"]
        MKDUP["samtools markdup"]
        BIDX["samtools index"]
        BAM["mapped/{sample}.bam<br/>+ .bai"]
    end

    FQ --> FASTP --> BWAMEM
    FASTP --> MULTIQC
    REF --> FAIDX
    REF --> BWAIDX
    BWAIDX --> BWAMEM
    BWAMEM --> MKDUP --> BIDX --> BAM
    BAM0 -.bam-input=True.-> BAM
    FAIDX --> BAM

    %% ====== Five callers in parallel ======
    subgraph CALL["CNV calling (5 tools, per-sample parallel)"]

        subgraph CK["CNVKit (split batch, flat reference, ploidy 1)"]
            CK_ACC["cnvkit access (cohort)"]
            CK_TGT["cnvkit target --split -a binSize (cohort)"]
            CK_ATG["cnvkit antitarget (cohort)"]
            CK_REF["cnvkit reference (flat, cohort)"]
            CK_TCOV["cnvkit coverage (target)"]
            CK_ACOV["cnvkit coverage (antitarget)"]
            CK_FIX["cnvkit fix --no-edge"]
            CK_SEG["cnvkit segment -m cbs -t 1e-6<br/>--drop-low-coverage"]
            CK_SMT["cnvkit segmetrics --ci --pi"]
            CK_CALL["cnvkit call -m clonal<br/>--purity 1 --ploidy 1 --filter ci"]
            CK_BED["res/cnvkit/{sample}.bed"]
            CK_ACC --> CK_TGT --> CK_ATG --> CK_REF
            CK_TGT --> CK_TCOV
            CK_ATG --> CK_ACOV
            CK_TCOV --> CK_FIX
            CK_ACOV --> CK_FIX
            CK_REF --> CK_FIX --> CK_SEG --> CK_SMT --> CK_CALL --> CK_BED
        end

        subgraph CP["CNVpytor (read-depth, haploid)"]
            CP_RUN["cnvpytor: -rd → -gc → -his → -partition → -call"]
            CP_BED["res/cnvpytor/{sample}.bed"]
            CP_RUN --> CP_BED
        end

        subgraph MO["cn.MOPS (haplocn.mops, cohort)"]
            MO_R["mopsCall.R (all BAMs at once)"]
            MO_BED["res/mops/{sample}.bed"]
            MO_R --> MO_BED
        end

        subgraph SM["Smoove (split-read SV)"]
            SM_CALL["smoove call (per-sample, single)"]
            SM_EX["bcftools query + DEL/DUP filter"]
            SM_BED["res/smoove/{sample}.bed"]
            SM_CALL --> SM_EX --> SM_BED
        end

        subgraph DL["Delly (SV mode, cohort genotyping)"]
            DL_CALL["delly call (per-sample)"]
            DL_MERGE["delly merge (cohort sites)"]
            DL_GENO["delly call -v merged (re-genotype per-sample)"]
            DL_GMERGE["bcftools merge"]
            DL_FILT["delly filter -f germline"]
            DL_UNC["bcftools view -s {sample}"]
            DL_EX["bcftools query + PASS + DEL/DUP"]
            DL_BED["res/delly/{sample}.bed"]
            DL_CALL --> DL_MERGE --> DL_GENO --> DL_GMERGE --> DL_FILT --> DL_UNC --> DL_EX --> DL_BED
        end
    end

    BAM --> CK_TCOV
    BAM --> CK_ACOV
    BAM --> CP_RUN
    BAM --> MO_R
    BAM --> SM_CALL
    BAM --> DL_CALL
    REF --> CK_ACC
    REF --> CP_RUN
    REF --> SM_CALL
    REF --> DL_CALL

    %% ====== Merge + Score + Filter + Annotate ======
    subgraph POST["Post-processing (per-sample)"]
        MERGE["mergeCNV.py<br/>priority: smoove → delly → cnvkit → cnvpytor → mops<br/>breakpoints locked to smoove/delly when present<br/>diploid→haploid CN conversion for smoove/delly"]
        BED2VCF["bed2vcf.py"]
        DUPHOLD["duphold (DHFC/DHBFC/DHFFC)"]
        DUPEX["bcftools query + scoreDuphold.py"]
        HF["hardFilter.py<br/>toolNum ≥ 2 · duphold > 0<br/>depth-only ≥ 600 bp"]
        ANNOT["annotate.py<br/>GS (low-complexity) + MS (low-mappability)"]
        FINAL["res/final/{sample}.bed"]
        MERGE --> BED2VCF --> DUPHOLD --> DUPEX --> HF --> ANNOT --> FINAL
    end

    CK_BED --> MERGE
    CP_BED --> MERGE
    MO_BED --> MERGE
    SM_BED --> MERGE
    DL_BED --> MERGE
    BAM --> DUPHOLD
    FAIDX --> BED2VCF
    REF --> DUPHOLD
    ANNO --> ANNOT

    %% ====== Visualization ======
    subgraph OUT["Cohort outputs"]
        PLOT["plotOverview.py<br/>res/report/{sample}.overview.png"]
    end

    FINAL --> PLOT

    %% Cross-sample population merge is a standalone post-workflow step:
    %%   filterCNV.py (optional refilter) -> mergeCNVPopulation.py -> cn/binary matrix

    %% ====== Styling ======
    classDef cohort fill:#fde8c4,stroke:#b97a00,color:#000;
    classDef persample fill:#d6ebff,stroke:#1f6feb,color:#000;
    classDef io fill:#eeeeee,stroke:#444,color:#000;
    classDef final fill:#d4f7d4,stroke:#1a7f37,color:#000;

    class FQ,BAM0,REF,ANNO io;
    class FAIDX,BWAIDX,MULTIQC,CK_ACC,CK_TGT,CK_ATG,CK_REF,MO_R,DL_MERGE,DL_GMERGE,DL_FILT cohort;
    class FASTP,BWAMEM,MKDUP,BIDX,BAM,CK_TCOV,CK_ACOV,CK_FIX,CK_SEG,CK_SMT,CK_CALL,CK_BED,CP_RUN,CP_BED,MO_BED,SM_CALL,SM_EX,SM_BED,DL_CALL,DL_GENO,DL_UNC,DL_EX,DL_BED,MERGE,BED2VCF,DUPHOLD,DUPEX,HF,ANNOT,PLOT persample;
    class FINAL final;
```

Legend: gray = inputs · orange = cohort-level (runs once across all samples) · blue = per-sample parallel · green = final outputs.

## Usage

```bash
# 1. Edit config.yaml: set genome, samples.tsv path, annotation BEDs, outdir
# 2. Edit samples.tsv: list sample names and fastq paths
# 3. Run
snakemake --profile profiles/default
```

The profile sets `--use-conda`, `--conda-prefix conda_envs`, `--latency-wait 60`, `-j 10`.

## Cross-sample population merge (standalone)

The workflow stops at per-sample `res/final/{sample}.bed`. Building the population
CN/binary matrices is a separate step so each species can use its own refilter and
overlap parameters. Two scripts, run after the workflow:

```bash
# 1. (optional) Refilter per-sample calls on duphold fold-change metrics.
#    DUP keeps metric > threshold; DEL keeps metric < threshold.
#    No threshold given for a type => that type passes through unchanged.
python scripts/filterCNV.py \
    --in-dir res/filtered --out-dir res/refilt_tier1 \
    --preset tier1                      # DUP all three >1.3, DEL all three <0.7

# fine-grained, e.g. constrain only the GC-corrected metric for DUPs:
python scripts/filterCNV.py \
    --in-dir res/filtered --out-dir res/refilt_dhbfc \
    --types dup --dup-dhbfc 1.5
# --mode all (default) requires every supplied threshold; --mode any needs one.
# Writes filter_log.tsv (parameters + per-sample kept counts) for reproducibility.

# 2. Cross-sample merge into population matrices.
python scripts/mergeCNVPopulation.py \
    --in-dir res/refilt_tier1 --out-dir res/population_tier1
# Overlap defaults match the old workflow (0.75 / 0.95 / 0.5, freq 0.8); override
# with --min-overlap/--max-overlap/--length-ratio/--freq-threshold.
```

Both read per-sample BEDs from `--in-dir` (default `res/filtered`), so you can merge
straight from hard-filtered calls or from a refiltered set. The merge is a greedy
front-to-back overlap merge, so region boundaries depend on input order; pass
`--samples` in a fixed order for reproducible boundaries across runs.

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
