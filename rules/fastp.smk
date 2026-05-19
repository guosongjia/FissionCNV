# =================================================================================================
#     Clean reads by fastp (paired-end only)
# =================================================================================================

def get_fastq(wildcards):
    """Get paired-end fastq files of given sample."""
    fastqs = config["global"]["samples"].loc[(wildcards.sample), ["fq1", "fq2"]]
    return {"r1": fastqs.fq1, "r2": fastqs.fq2}

def unpack_fastq_files(wildcards):
    return list(get_fastq(wildcards).values())

rule clean_reads_pe:
    input:
        unpack_fastq_files,
    output:
        trimmed = (
            ["cleaned/{sample}_1.fq.gz", "cleaned/{sample}_2.fq.gz"]
            if config['settings']['keep-intermediate']['fastp']
            else temp(["cleaned/{sample}_1.fq.gz", "cleaned/{sample}_2.fq.gz"])
        ),
        html = "cleaned/{sample}-pe-fastp.html",
        json = "cleaned/{sample}-pe-fastp.json",
    log:
        "logs/fastp/{sample}.log"
    benchmark:
        "benchmarks/fastp/{sample}.bench"
    threads:
        config["params"]["fastp"]["threads"]
    conda:
        "../envs/preprocess.yaml"
    shell:
        "(fastp --thread {threads} --in1 {input[0]} --in2 {input[1]} "
        "--out1 {output.trimmed[0]} --out2 {output.trimmed[1]} "
        "--html {output.html} --json {output.json}) > {log} 2>&1"

rule multiqc_report:
    input:
        json = expand("cleaned/{sample}-pe-fastp.json", sample=config["global"]["sample-names"]),
    output:
        "cleaned/multiqc-report.html",
    conda:
        "../envs/preprocess.yaml"
    shell:
        "multiqc --force -d cleaned/ -n multiqc-report -o cleaned/ -q"

localrules: all_fastp

rule all_fastp:
    input:
        cleaned_reads = expand("cleaned/{sample}_{pair}.fq.gz", pair=[1, 2],
            sample=config["global"]["sample-names"]),
        report = "cleaned/multiqc-report.html",

def get_cleaned_reads(wildcards):
    if config['params']['skip-fastp']:
        fastqs = config["global"]["samples"].loc[(wildcards.sample), ["fq1", "fq2"]]
        return [fastqs.fq1, fastqs.fq2]
    return expand("cleaned/{sample}_{pair}.fq.gz", pair=[1, 2], sample=wildcards.sample)
