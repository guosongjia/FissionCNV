# =================================================================================================
#     Visualization: per-sample whole-genome CNV overview
# =================================================================================================

rule plot_overview:
    input:
        "res/final/{sample}.bed",
    output:
        "res/report/{sample}.overview.png",
    params:
        fai = config['data']['genome'] + ".fai",
    conda:
        "../envs/postprocess.yaml"
    script:
        "../scripts/plotOverview.py"
