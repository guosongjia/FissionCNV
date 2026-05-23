# =================================================================================================
#     Use duphold to measure read depth of CNV adjacent regions
# =================================================================================================

rule convert_bed2vcf:
    input:
        bed = "res/merge/{sample}.merged.bed",
        fai = config['data']['genome'] + ".fai",
    output:
        temp("temp/duphold/{sample}.vcf"),
    params:
        absPath = config['params']['absPath']
    shell:
        "python {params.absPath}/scripts/bed2vcf.py {input.bed} {input.fai} {output}"

rule duphold_score:
    input:
        "mapped/{sample}.bam.bai",
        bam = "mapped/{sample}.bam",
        vcf = rules.convert_bed2vcf.output,
    output:
        temp("temp/duphold/{sample}.duphold.vcf"),
    params:
        genome = config['data']['genome'],
    threads: 8
    log:
        "logs/duphold/{sample}.duphold.log"
    benchmark:
        "benchmarks/duphold/{sample}.duphold.bench"
    conda:
        "../envs/postprocess.yaml"
    shell:
        "duphold -t {threads} -v {input.vcf} -b {input.bam} -f {params.genome} -o {output}"

rule duphold_extract:
    input:
        rules.duphold_score.output,
    output:
        temp("temp/duphold/{sample}.duphold.bed"),
    conda:
        "../envs/postprocess.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%CN\t%AS\t%DHFC\t%DHBFC\t%DHFFC]\t%INFO/TNa\t%INFO/TN\t%INFO/GCF\t%INFO/SAMPLE\n' {input} > {output}"

rule duphold_convert:
    input:
        rules.duphold_extract.output,
    output:
        "res/duphold/{sample}.bed",
    script:
        "../scripts/scoreDuphold.py"
