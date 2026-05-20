# =================================================================================================
#     CNV calling by Smoove (split-read based)
# =================================================================================================

rule smoove_call:
    input:
        "mapped/{sample}.bam.bai",
        bam = "mapped/{sample}.bam",
    output:
        "res/smoove/{sample}.raw.vcf.gz",
    params:
        outdir = "res/smoove/",
        ref = config['data']['genome'],
    log:
        "logs/smoove/{sample}.call.log"
    benchmark:
        "benchmarks/smoove/{sample}.bench"
    conda:
        "../envs/smoove.yaml"
    shell:
        "(smoove call --outdir {params.outdir} "
        "--name {wildcards.sample} --fasta {params.ref} -p 1 {input.bam} && "
        "mv {params.outdir}{wildcards.sample}-smoove.vcf.gz {output}) > {log} 2>&1"

rule smoove_extract:
    input:
        rules.smoove_call.output,
    output:
        "temp/smoove/{sample}.bed",
    conda:
        "../envs/smoove.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL\n' {input} | "
        "egrep 'DUP|DEL' | awk -v OFS='\t' '$4==\"DEL\" {{print $1,$2,$3,1,$5}} "
        "$4==\"DUP\" {{print $1,$2,$3,3,$5}}' > {output}"

rule smoove_convert:
    input:
        rules.smoove_extract.output,
    output:
        "res/smoove/{sample}.bed",
    script:
        "../scripts/smooveConvert.py"

localrules: all_smoove

rule all_smoove:
    input:
        expand("res/smoove/{sample}.bed", sample=config['global']['sample-names'])
