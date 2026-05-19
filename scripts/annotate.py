#! /usr/bin/env python
# Annotate CNVs with GS (low-complexity) and MS (low-mappability) scores.
# These are annotation columns only, not used for hard filtering.

import sys
sys.path.insert(0, snakemake.params['absPath'] + '/scripts')
from utils import readCNVFile, calculateOverlapProp4Region


def calculateOverlapScore(target_cnv, cnv_list, reciprocal_prop=0.3):
    """Calculate overlap score: 100 = no overlap with bad regions."""
    accumLen = 0
    i = 0
    for cnv in cnv_list:
        overlap, prop1, prop2 = calculateOverlapProp4Region(target_cnv, cnv)
        if min(prop1, prop2) > reciprocal_prop:
            i += 1
            accumLen += overlap
    cnv_len = int(target_cnv[2]) - int(target_cnv[1])
    if cnv_len == 0:
        return 100
    accumProp = round(accumLen * 100 / cnv_len)
    score = 100 - accumProp - i * 2
    return score


inputFile = snakemake.input[0]
outputFile = snakemake.output[0]
lowComplexFile = snakemake.params['low_complexity']
lowMapFile = snakemake.params['low_mappable']

# Load annotation BED files (if provided)
bad_list = readCNVFile(lowComplexFile, tool='Bad') if lowComplexFile else []
lowMap_list = readCNVFile(lowMapFile, tool='Bad') if lowMapFile else []

with open(inputFile, 'r') as f, open(outputFile, 'w') as g:
    header = f.readline().strip()
    print(header + '\tGS\tMS', file=g)
    for line in f:
        cnv = line.strip().split('\t')[:3]
        gs = calculateOverlapScore(cnv, bad_list) if bad_list else 100
        ms = calculateOverlapScore(cnv, lowMap_list) if lowMap_list else 100
        print(line.strip(), gs, ms, sep='\t', file=g)
