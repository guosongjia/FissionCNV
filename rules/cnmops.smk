# =================================================================================================
#     CNV calling by cn.MOPS (haplocn.mops for haploid)
# =================================================================================================

# cn.MOPS requires at least 6 samples to run as it calls CNVs based on a mixed Poisson model.
# Using haplocn.mops() for haploid fission yeast genomes.
rule mops_call:
    input:
        get_sample_bai(config['global']['sample-names']),
        bam = get_sample_bam(config['global']['sample-names']),
    output:
        bed = expand("res/mops/{sample}.raw.bed", sample=config['global']['sample-names']),
    params:
        resDir = "res/mops/",
        binSize = config['params']['binSize'],
        absPath = config['params']['absPath'],
    threads:
        config['params']['mops']['threads']
    log:
        "logs/mops/call.log"
    benchmark:
        "benchmarks/mops/call.bench"
    conda:
        "../envs/cnmops.yaml"
    shell:
        "Rscript {params.absPath}/scripts/mopsCall.R {params.resDir} {params.binSize} "
        "{threads} {input.bam} > {log} 2>&1"

rule mops_convert:
    input:
        "res/mops/{sample}.raw.bed",
    output:
        "res/mops/{sample}.bed",
    params:
        absPath = config['params']['absPath'],
    script:
        "../scripts/mopsConvert.py"

localrules: all_mops

rule all_mops:
    input:
        expand("res/mops/{sample}.bed", sample=config['global']['sample-names']),
