# =================================================================================================
#     Pre-processing
# =================================================================================================

# Build index for reference genome
rule samtools_faidx:
    input:
        config["data"]["genome"]
    output:
        config["data"]["genome"] + ".fai"
    conda:
        "../envs/preprocess.yaml"
    shell:
        "samtools faidx {input}"

# Build bwa index for reference genome
if not config['params']['bwamem']['index']:
    rule bwa_index:
        input:
            config["data"]["genome"]
        output:
            [config['data']['genome'] + ext for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
        log:
            "logs/index/bwa_index.log"
        conda:
            "../envs/preprocess.yaml"
        shell:
            "bwa index {input} > {log} 2>&1"

localrules: all_prep

rule all_prep:
    input:
        ref = config["data"]["genome"],
        ref_idx = [config['data']['genome'] + ext for ext in [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]],
