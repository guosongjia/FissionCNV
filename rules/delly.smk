# =================================================================================================
#     CNV calling by Delly (SV mode for non-human species)
# =================================================================================================

# Delly calls structural variants, then merges, genotypes, and filters across samples.
rule delly_call_sv:
    input:
        "mapped/{sample}.bam.bai",
        bam = "mapped/{sample}.bam",
    output:
        "temp/delly/call/{sample}.sv.bcf",
    params:
        ref = config['data']['genome'],
    log:
        "logs/delly/{sample}.callsv.log"
    benchmark:
        "benchmarks/delly/{sample}.callsv.bench"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly call -g {params.ref} -o {output} {input.bam} > {log} 2>&1"

rule delly_merge:
    input:
        expand("temp/delly/call/{sample}.sv.bcf", sample=config['global']['sample-names']),
    output:
        "temp/delly/mergedSites.sv.bcf",
    log:
        "logs/delly/mergeSV.log"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly merge -o {output} {input} > {log} 2>&1"

rule delly_genotype:
    input:
        bam = "mapped/{sample}.bam",
        merged = "temp/delly/mergedSites.sv.bcf",
    output:
        "temp/delly/genotype/{sample}.geno.sv.bcf",
    params:
        ref = config['data']['genome'],
    log:
        "logs/delly/{sample}.genotype.log"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly call -g {params.ref} -v {input.merged} -o {output} {input.bam} > {log} 2>&1"

rule delly_genotype_merge:
    input:
        expand("temp/delly/genotype/{sample}.geno.sv.bcf", sample=config['global']['sample-names']),
    output:
        "temp/delly/genotype.merged.bcf",
    threads: 10
    log:
        "logs/delly/genotype.merge.log"
    conda:
        "../envs/delly.yaml"
    shell:
        "bcftools merge --threads {threads} -m id -O b -o {output} {input} > {log} 2>&1; "
        "bcftools index {output}"

rule delly_filter:
    input:
        "temp/delly/genotype.merged.bcf",
    output:
        "temp/delly/germline.bcf"
    log:
        "logs/delly/filter.log"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly filter -f germline -o {output} {input} > {log} 2>&1"

rule delly_uncompress:
    input:
        rules.delly_filter.output,
    output:
        "res/delly/{sample}.raw.bcf",
    conda:
        "../envs/delly.yaml"
    shell:
        "bcftools view -s {wildcards.sample} -O b -o {output} {input}"

rule delly_extract:
    input:
        rules.delly_uncompress.output,
    output:
        "temp/delly/extract/{sample}.bed",
    conda:
        "../envs/delly.yaml"
    shell:
        "bcftools query -f '[%FT]\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%RDCN]\n' {input} | "
        "grep 'PASS' | egrep 'DEL|DUP' | cut -f 2- | awk '$4 == \"DEL\" && $5 < 2 {{print$0}} "
        "$4 == \"DUP\" && $5 > 2 {{print$0}}' | cut -f 1,2,3,5 > {output}"

rule delly_convert:
    input:
        rules.delly_extract.output,
    output:
        "res/delly/{sample}.bed",
    script:
        "../scripts/smooveConvert.py"

localrules: all_delly

rule all_delly:
    input:
        expand("res/delly/{sample}.bed", sample=config['global']['sample-names'])
