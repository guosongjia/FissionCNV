#! /usr/bin/env python
# Hard filter CNVs based on:
# 1. toolNum >= min_tools
# 2. dupholdScore > min_duphold_score
# 3. If only depth-based tools support: length >= depth_only_min_length

inputFile = snakemake.input[0]
outputFile = snakemake.output[0]
min_tools = snakemake.params['min_tools']
depth_only_min_length = snakemake.params['depth_only_min_length']

SPLITREAD_TOOLS = {'smoove', 'delly'}

with open(inputFile, 'r') as f, open(outputFile, 'w') as g:
    header = f.readline().strip()
    print(header, file=g)
    for line in f:
        x = line.strip().split('\t')
        start, end = int(x[1]), int(x[2])
        toolNum = int(x[11])
        dupholdScore = int(x[6])
        toolName = x[10]
        length = end - start

        # Hard filter 1: minimum tool support
        if toolNum < min_tools:
            continue
        # Hard filter 2: duphold score > 0
        if dupholdScore <= 0:
            continue
        # Hard filter 3: depth-only calls must be >= 600bp
        tools_set = set(toolName.split(','))
        has_splitread = bool(tools_set & SPLITREAD_TOOLS)
        if not has_splitread and length < depth_only_min_length:
            continue

        print(line.strip(), file=g)
