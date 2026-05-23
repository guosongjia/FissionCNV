# =================================================================================================
#     CNV calling by CNVKit (flat reference, haploid)
# =================================================================================================
# Original `cnvkit batch` is split into cohort-level reference prep (4 rules)
# + per-sample pipeline (4 rules) so Snakemake can parallelize per-sample work.
# Functionally equivalent to: cnvkit batch -n -m wgs --target-avg-size {binSize}
#                             -f genome.fa --drop-low-coverage

# -----------------------------------------------------------------------------
# Cohort-level reference prep (runs once per genome)
# -----------------------------------------------------------------------------
rule cnvkit_access:
    input:
        config['data']['genome'],
    output:
        temp("temp/cnvkit/access.bed"),
    log:
        "logs/cnvkit/access.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py access {input} -o {output} > {log} 2>&1"

rule cnvkit_target:
    input:
        rules.cnvkit_access.output,
    output:
        temp("temp/cnvkit/target.bed"),
    params:
        binSize = config['params']['binSize'],
    log:
        "logs/cnvkit/target.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py target {input} --split -a {params.binSize} -o {output} > {log} 2>&1"

rule cnvkit_antitarget:
    input:
        access = rules.cnvkit_access.output,
        target = rules.cnvkit_target.output,
    output:
        temp("temp/cnvkit/antitarget.bed"),
    log:
        "logs/cnvkit/antitarget.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py antitarget {input.target} --access {input.access} -o {output} > {log} 2>&1"

rule cnvkit_reference:
    input:
        target     = rules.cnvkit_target.output,
        antitarget = rules.cnvkit_antitarget.output,
        fasta      = config['data']['genome'],
    output:
        temp("temp/cnvkit/myFlatReference.cnn"),
    log:
        "logs/cnvkit/reference.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py reference -t {input.target} -a {input.antitarget} "
        "-f {input.fasta} -o {output} > {log} 2>&1"

# -----------------------------------------------------------------------------
# Per-sample pipeline (runs in parallel up to Snakemake `cores` budget)
# -----------------------------------------------------------------------------
rule cnvkit_target_coverage:
    input:
        bai    = "mapped/{sample}.bam.bai",
        bam    = "mapped/{sample}.bam",
        target = rules.cnvkit_target.output,
    output:
        temp("temp/cnvkit/{sample}.targetcoverage.cnn"),
    threads: 1
    log:
        "logs/cnvkit/coverage/{sample}.target.log",
    benchmark:
        "benchmarks/cnvkit/{sample}.target_coverage.bench"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py coverage {input.bam} {input.target} -o {output} > {log} 2>&1"

rule cnvkit_antitarget_coverage:
    input:
        bai        = "mapped/{sample}.bam.bai",
        bam        = "mapped/{sample}.bam",
        antitarget = rules.cnvkit_antitarget.output,
    output:
        temp("temp/cnvkit/{sample}.antitargetcoverage.cnn"),
    threads: 1
    log:
        "logs/cnvkit/coverage/{sample}.antitarget.log",
    benchmark:
        "benchmarks/cnvkit/{sample}.antitarget_coverage.bench"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py coverage {input.bam} {input.antitarget} -o {output} > {log} 2>&1"

rule cnvkit_fix:
    input:
        target     = rules.cnvkit_target_coverage.output,
        antitarget = rules.cnvkit_antitarget_coverage.output,
        reference  = rules.cnvkit_reference.output,
    output:
        temp("temp/cnvkit/{sample}.cnr"),
    log:
        "logs/cnvkit/fix/{sample}.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py fix {input.target} {input.antitarget} {input.reference} "
        "--no-edge -o {output} > {log} 2>&1"

rule cnvkit_segment:
    input:
        rules.cnvkit_fix.output,
    output:
        temp("temp/cnvkit/{sample}.cns"),
    threads: 1
    log:
        "logs/cnvkit/segment/{sample}.log",
    benchmark:
        "benchmarks/cnvkit/{sample}.segment.bench"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py segment {input} -m cbs -t 1e-6 --drop-low-coverage "
        "-o {output} > {log} 2>&1"

# -----------------------------------------------------------------------------
# Segmetrics (CI for filter), call (haploid clonal), and BED conversion
# -----------------------------------------------------------------------------
rule cnvkit_segmetric:
    input:
        cns = rules.cnvkit_segment.output,
        cnr = rules.cnvkit_fix.output,
    output:
        temp("temp/cnvkit/segmetrics/{sample}.cns"),
    params:
        "--ci --pi"
    log:
        "logs/cnvkit/segmetric/{sample}.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py segmetrics -s {input.cns} {input.cnr} {params} -o {output} > {log} 2>&1"

rule cnvkit_call:
    input:
        rules.cnvkit_segmetric.output,
    output:
        "res/cnvkit/{sample}.raw.cns",
    params:
        "-m clonal --purity 1 --ploidy 1 --filter ci"
    log:
        "logs/cnvkit/call/{sample}.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py call {input} {params} -o {output} > {log} 2>&1"

rule cnvkit_convert:
    input:
        rules.cnvkit_call.output,
    output:
        "res/cnvkit/{sample}.bed",
    params:
        absPath = config['params']['absPath'],
    script:
        "../scripts/cnvkitConvert.py"

localrules: all_cnvkit

rule all_cnvkit:
    input:
        expand("res/cnvkit/{sample}.bed", sample=config['global']['sample-names']),
