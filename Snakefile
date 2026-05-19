#! /usr/bin/env snakemake

configfile: "config.yaml"

workdir: config['data']['outdir']

include: "rules/common.smk"

# Preprocessing: QC and alignment (skip if BAM input)
if not config['params']['bam-input']:
    include: "rules/pre-processing.smk"
    include: "rules/fastp.smk"
    include: "rules/bwamem.smk"

# CNV calling by 5 methods
include: "rules/cnvkit.smk"
include: "rules/cnvpytor.smk"
include: "rules/cnmops.smk"
include: "rules/smoove.smk"
include: "rules/delly.smk"

# Merge, filter, annotate, and generate population matrix
include: "rules/merge.smk"

# Visualization
include: "rules/report.smk"

rule all:
    input:
        expand("res/final/{sample}.bed", sample=config['global']['sample-names']),
        "res/population/cn_matrix.tsv",
        "res/population/binary_matrix.tsv",
        expand("res/report/{sample}.overview.png", sample=config['global']['sample-names']),
