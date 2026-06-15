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

# Merge tool calls, filter, and annotate (cross-sample population merge is a
# standalone step: scripts/filterCNV.py + scripts/mergeCNVPopulation.py)
include: "rules/merge.smk"

# Breakpoint analysis (optional)
if config['params'].get('breakpoint-analysis', False):
    include: "rules/breakpoint.smk"

# Visualization
include: "rules/report.smk"

rule all:
    input:
        expand("res/final/{sample}.bed", sample=config['global']['sample-names']),
        expand("res/report/{sample}.overview.png", sample=config['global']['sample-names']),
