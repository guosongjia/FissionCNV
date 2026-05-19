#! /usr/bin/env python

# Convert cnvpytor output into bed format CNV list for haploid genome.

inputfile = snakemake.input[0]
outputfile = snakemake.output[0]

with open(inputfile, 'r') as f:
    with open(outputfile, 'w') as g:
        print('chromosome', 'start', 'end', 'cn', 'log2', 'evalue|pN|dG', sep='\t', file=g)
        for line in f:
            line = line.strip().split('\t')
            coord = line[1].split(':')
            chrom = coord[0]
            start = coord[1].split('-')[0]
            end = coord[1].split('-')[1]
            depth = float(line[3])
            cn = round(depth)  # haploid: CN = normalized depth ratio
            evalue1 = float(line[4])
            evalue2 = float(line[5])
            evalue = max(evalue1, evalue2)
            q0 = float(line[8])
            pN = float(line[9])
            dG = int(line[10])
            if cn != 1 and evalue < 0.00001 and q0 < 0.5 and pN < 0.5 and dG > 1000:
                print(chrom, start, end, cn, depth, '|'.join(str(x) for x in [evalue, pN, dG]),
                sep='\t', file=g)
