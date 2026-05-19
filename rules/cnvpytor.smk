# =================================================================================================
#     CNV calling by cnvpytor (haploid)
# =================================================================================================

# CNVpytor calls CNVs sample by sample, does not require control sample.
# For fission yeast, call all chromosomes (no species-specific filtering).
rule cnvpytor_call:
    input:
        "mapped/{sample}.bam.bai",
        bam = "mapped/{sample}.bam",
    output:
        pytor = "temp/cnvpytor/{sample}.pytor",
        call = "temp/cnvpytor/{sample}.call",
    params:
        binSize = config['params']['binSize'],
        ref = config['data']['genome'],
    threads:
        config['params']['cnvpytor']['threads']
    log:
        "logs/cnvpytor/{sample}.call.log"
    benchmark:
        "benchmarks/cnvpytor/{sample}.call.bench"
    conda:
        "../envs/cnvpytor.yaml"
    shell:
        "(cnvpytor -root {output.pytor} -j {threads} -rd {input.bam}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -gc {params.ref} -make_gc_file; \n"
        "cnvpytor -root {output.pytor} -j {threads} -his {params.binSize}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -partition {params.binSize}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -call {params.binSize} > {output.call}) > {log} 2>&1"

rule cnvpytor_convert:
    input:
        "temp/cnvpytor/{sample}.call",
    output:
        "res/cnvpytor/{sample}.bed",
    script:
        "../scripts/cnvpytorConvert.py"

localrules: all_cnvpytor

rule all_cnvpytor:
    input:
        expand("res/cnvpytor/{sample}.bed", sample=config['global']['sample-names']),
