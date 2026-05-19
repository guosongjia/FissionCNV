# =================================================================================================
#     Mapping with BWA MEM
# =================================================================================================

def get_bwa_mem_rg(wildcards):
    rg = ["ID:" + wildcards.sample, "SM:" + wildcards.sample, "PL:" + config['params']['bwamem']['platform']]
    tmp_rg_tags = "\\t".join(rg)
    rg_tags = "-R '@RG\\t" + tmp_rg_tags + "' "
    return rg_tags

rule map_reads:
    input:
        reads = get_cleaned_reads,
        idx = [config['data']['genome'] + ext for ext in [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]],
        ref = config['data']['genome'],
    output:
        ("mapped/{sample}.raw.bam" if config['settings']['keep-intermediate']['bwamem']
            else temp("mapped/{sample}.raw.bam"))
    params:
        read_group = get_bwa_mem_rg,
    threads:
        config['params']['bwamem']['threads']
    log:
        "logs/bwamem/{sample}.log"
    benchmark:
        "benchmarks/bwamem/{sample}.bench"
    conda:
        "../envs/preprocess.yaml"
    shell:
        "(bwa mem -M {params.read_group} -t {threads} {input.ref} {input.reads} | "
        "samtools fixmate -m -@ 4 - - | "
        "samtools sort -@ 10 -o {output}) >{log} 2>&1"

# Mark duplicates (mark only, do not remove)
rule markDuplicates:
    input:
        rules.map_reads.output,
    output:
        bam = "mapped/{sample}.bam",
    threads: 6
    log:
        "logs/markdup/{sample}.log"
    conda:
        "../envs/preprocess.yaml"
    shell:
        "samtools markdup -@ {threads} {input} {output.bam} > {log} 2>&1"

rule samtools_index:
    input:
        "mapped/{sample}.bam",
    output:
        "mapped/{sample}.bam.bai",
    log:
        "logs/samtools/{sample}.index.log",
    conda:
        "../envs/preprocess.yaml"
    shell:
        "samtools index {input} > {log} 2>&1"

localrules: all_bwamem

rule all_bwamem:
    input:
        expand("mapped/{sample}.bam", sample=config['global']['sample-names']),
        expand("mapped/{sample}.bam.bai", sample=config['global']['sample-names'])
