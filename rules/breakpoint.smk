# =================================================================================================
#     Breakpoint annotation: repeat context + BAM evidence + microhomology (precise tools only)
# =================================================================================================

rule breakpoint_annotate:
    input:
        annotated_bed = "res/annotated/{sample}.bed",
        bam = "mapped/{sample}.bam",
        bai = "mapped/{sample}.bam.bai",
        smoove_vcf = "res/smoove/{sample}.raw.vcf.gz",
        delly_bcf = "res/delly/{sample}.raw.bcf",
    output:
        "res/final/{sample}.bed",
    params:
        absPath = config['params']['absPath'],
        genome = config['data']['genome'],
        low_complexity = config['data'].get('low-complexity', ''),
        trna_gff = config['data'].get('breakpoint-annotations', {}).get('trna_gff', ''),
        rrna_gff = config['data'].get('breakpoint-annotations', {}).get('rrna_gff', ''),
        ltr_bed = config['data'].get('breakpoint-annotations', {}).get('ltr_bed', ''),
        centromeric_bed = config['data'].get('breakpoint-annotations', {}).get('centromeric_bed', ''),
        subtelomeric_bed = config['data'].get('breakpoint-annotations', {}).get('subtelomeric_bed', ''),
        kmds_bed = config['data'].get('breakpoint-annotations', {}).get('kmds_bed', ''),
    conda:
        "../envs/breakpoint.yaml"
    script:
        "../scripts/breakpointAnnotate.py"
