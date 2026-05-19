#! /usr/bin/env python

# Merge CNV calls from 5 tools for haploid fission yeast.
# Priority order: smoove -> delly -> cnvkit -> cnvpytor -> mops
# Smoove/Delly output diploid CN values; convert to haploid before merging.

import sys
sys.path.insert(0, snakemake.params['absPath'] + '/scripts')
from utils import readCNVFile, mergeCNVFromTools, convert_diploid_cn_to_haploid

sample = snakemake.wildcards.sample
absPath = snakemake.params['absPath']
outputFile = snakemake.output[0]
min_overlap = snakemake.params['min_overlap']
max_overlap = snakemake.params['max_overlap']

# Read CNVs from each tool in priority order
tools_order = ['smoove', 'delly', 'cnvkit', 'cnvpytor', 'mops']
tool_files = {
    'smoove': f'res/smoove/{sample}.bed',
    'delly': f'res/delly/{sample}.bed',
    'cnvkit': f'res/cnvkit/{sample}.bed',
    'cnvpytor': f'res/cnvpytor/{sample}.bed',
    'mops': f'res/mops/{sample}.bed',
}

# Tools that output diploid CN values
diploid_tools = {'smoove', 'delly'}

cnvList = []
tool_type_map = {
    'smoove': 'Smoove',
    'delly': 'Delly',
    'cnvkit': 'MOPS',
    'cnvpytor': 'CNVpytor',
    'mops': 'MOPS',
}
for tool in tools_order:
    tool_type = tool_type_map[tool]
    cnvs = readCNVFile(tool_files[tool], tool=tool_type)
    for cnv in cnvs:
        if tool in diploid_tools:
            cnv[3] = convert_diploid_cn_to_haploid(cnv[3])
        if cnv[3] != 1:  # skip normal CN after conversion
            cnv.append(tool)
            cnvList.append(cnv)

# Merge CNVs from all tools
mergedCnvs = mergeCNVFromTools(cnvList, min_threshold=min_overlap, max_threshold=max_overlap)

# Output: chrom, start, end, cn, tools, toolNum, accumScore
with open(outputFile, 'w') as f:
    print('chromosome', 'start', 'end', 'cn', 'tools', 'toolNum', 'accumScore', 'sample', sep='\t', file=f)
    for cnv in mergedCnvs:
        print(cnv[0], cnv[1], cnv[2], cnv[3], cnv[4], cnv[5], cnv[6], sample, sep='\t', file=f)
