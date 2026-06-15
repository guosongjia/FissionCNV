# =================================================================================================
#     Merge CNVs, score with duphold, filter, and annotate
# =================================================================================================

# Include duphold rules
include: "duphold.smk"

# Step 1: Merge CNV calls from 5 tools
rule merge_CNVCall:
    input:
        smoove = "res/smoove/{sample}.bed",
        delly = "res/delly/{sample}.bed",
        cnvkit = "res/cnvkit/{sample}.bed",
        cnvpytor = "res/cnvpytor/{sample}.bed",
        mops = "res/mops/{sample}.bed",
    output:
        "res/merge/{sample}.merged.bed",
    params:
        absPath = config['params']['absPath'],
        min_overlap = config['params']['merge']['min-reciprocal-overlap'],
        max_overlap = config['params']['merge']['max-overlap'],
    script:
        "../scripts/mergeCNV.py"

# Step 2: Duphold scoring (defined in duphold.smk)
# convert_bed2vcf -> duphold_score -> duphold_extract -> duphold_convert
# Output: res/duphold/{sample}.bed

# Step 3: Hard filtering
rule hard_filter:
    input:
        "res/duphold/{sample}.bed",
    output:
        "res/filtered/{sample}.bed",
    params:
        min_tools = config['params']['filter']['min-tools'],
        depth_only_min_length = config['params']['filter']['depth-only-min-length'],
    script:
        "../scripts/hardFilter.py"

# Step 4: Annotate with GS and MS scores
rule annotate:
    input:
        "res/filtered/{sample}.bed",
    output:
        "res/annotated/{sample}.bed",
    params:
        absPath = config['params']['absPath'],
        low_complexity = config['data']['low-complexity'],
        low_mappable = config['data']['low-mappable'],
    script:
        "../scripts/annotate.py"

# Step 5 (no-breakpoint path): passthrough to res/final when breakpoint-analysis is disabled
if not config['params'].get('breakpoint-analysis', False):
    rule final_passthrough:
        input:
            "res/annotated/{sample}.bed",
        output:
            "res/final/{sample}.bed",
        shell:
            "cp {input} {output}"

# Cross-sample merge (population CN / binary matrices) is intentionally NOT part
# of the workflow: different species often need different per-sample refilter and
# overlap parameters before merging. Run it standalone after the workflow finishes:
#   python scripts/filterCNV.py        --in-dir res/filtered --out-dir res/refilt ...
#   python scripts/mergeCNVPopulation.py --in-dir res/refilt   --out-dir res/population ...
