#! /usr/bin/env python

# Includes all functions used in FissionCNV
import os, subprocess
import pandas as pd

# Normal copy number for haploid
NORMAL_CN = 1


def testSameCNVType(cn1, cn2):
    """Test if two CNVs are the same type (both DEL or both DUP)
    >>> testSameCNVType(0, 0)
    True
    >>> testSameCNVType(2, 3)
    True
    >>> testSameCNVType(0, 2)
    False
    """
    if (cn1 > NORMAL_CN and cn2 > NORMAL_CN) or (cn1 < NORMAL_CN and cn2 < NORMAL_CN):
        return True
    return False


def readCNVFile(infile, tool):
    """Read file content of different sources into a CNV list"""
    cnvList = []
    with open(infile) as f:
        for line in f:
            if line.startswith('chrom'):
                continue
            x = line.strip().split('\t')
            if tool == 'CNVKit':
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[5])
            elif tool == 'CNVpytor':
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
            elif tool in ['MOPS', 'merge']:
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
            elif tool in ['Smoove', 'Delly']:
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), int(x[3])
            elif tool in ['Bad', 'Normal']:
                chrom, start, end, cn = x[0], int(x[1]), int(x[2]), 0
            else:
                raise ValueError('No tool is indicated, please indicate one!')
            if cn == NORMAL_CN:
                continue
            cnvList.append([chrom, start, end, cn])
    return cnvList


def convert_diploid_cn_to_haploid(cn):
    """Convert diploid CN (from smoove/delly) to haploid CN.
    >>> convert_diploid_cn_to_haploid(1)
    0
    >>> convert_diploid_cn_to_haploid(3)
    2
    >>> convert_diploid_cn_to_haploid(4)
    2
    """
    return round(cn / 2)


def testOverlapped(cnv1, cnv2):
    """ Test whether two CNVs are overlapped
    >>> testOverlapped(['chr1',1,10], ['chr1',10,20])
    False
    >>> testOverlapped(['chr1',5,15], ['chr1',10,20])
    True
    >>> testOverlapped(['chr2',10,20], ['chr1',10,20])
    False
    """
    c1, s1, e1 = cnv1[0], int(cnv1[1]), int(cnv1[2])
    c2, s2, e2 = cnv2[0], int(cnv2[1]), int(cnv2[2])
    if c1 == c2:
        if s2 < s1 < e2 or s2 < e1 < e2 or s1 < s2 < e2 < e1:
            return True
    return False


def calculateOverlapProp(cnv1, cnv2):
    """ Calculate the overlapped proportion between two CNVs
    >>> calculateOverlapProp(['chr1',100,1000,0], ['chr1',500,1000,0])
    (500, 0.56, 1.0)
    """
    c1, s1, e1, cn1 = cnv1[0], int(cnv1[1]), int(cnv1[2]), int(cnv1[3])
    c2, s2, e2, cn2 = cnv2[0], int(cnv2[1]), int(cnv2[2]), int(cnv2[3])
    cnvLen1, cnvLen2 = e1 - s1, e2 - s2
    overlap = 0
    assert cnvLen1 > 0 and cnvLen2 > 0
    if c1 == c2:
        if testSameCNVType(cn1, cn2):
            if s2 <= s1 < e1 <= e2:
                overlap = e1 - s1
            elif s2 <= s1 <= e2 < e1:
                overlap = e2 - s1
            elif s1 <= s2 < e2 <= e1:
                overlap = e2 - s2
            elif s1 < s2 <= e1 <= e2:
                overlap = e1 - s2
    prop1 = round((overlap / cnvLen1), 2)
    prop2 = round((overlap / cnvLen2), 2)
    return overlap, prop1, prop2


def calculateOverlapProp4Region(cnv1, cnv2):
    """Calculate the overlapped proportion between CNV and genomic region"""
    c1, s1, e1 = cnv1[0], int(cnv1[1]), int(cnv1[2])
    c2, s2, e2 = cnv2[0], int(cnv2[1]), int(cnv2[2])
    cnvLen1, cnvLen2 = e1 - s1, e2 - s2
    overlap = 0
    assert cnvLen1 > 0 and cnvLen2 > 0
    if c1 == c2:
        if s2 <= s1 < e1 <= e2:
            overlap = e1 - s1
        elif s2 <= s1 <= e2 < e1:
            overlap = e2 - s1
        elif s1 <= s2 < e2 <= e1:
            overlap = e2 - s2
        elif s1 < s2 <= e1 <= e2:
            overlap = e1 - s2
    prop1 = round((overlap / cnvLen1), 2)
    prop2 = round((overlap / cnvLen2), 2)
    return overlap, prop1, prop2


def mergeConsecutiveSegments(cnvList, shift=0):
    """ Merge consecutive segments and average copy number
    >>> mergeConsecutiveSegments([['chr1',1,10,0], ['chr1',10,20,0], ['chr2',1,5,2], ['chr2',5,10,3]], shift=0)
    [['chr1', 1, 20, 0], ['chr2', 1, 10, 2]]
    """
    new_cnvList = []
    if len(cnvList) == 0:
        return new_cnvList
    tmpChrom, tmpStart, tmpEnd, tmpCN = cnvList[0]
    for i in range(1, len(cnvList)):
        chrom, start, end, cn = cnvList[i]
        if chrom == tmpChrom and testSameCNVType(tmpCN, cn) and start == tmpEnd + shift:
            tmpEnd = end
            tmpCN = round((tmpCN + cn) / 2)
        else:
            new_cnvList.append([tmpChrom, tmpStart, tmpEnd, tmpCN])
            tmpChrom, tmpStart, tmpEnd, tmpCN = chrom, start, end, cn
    new_cnvList.append([tmpChrom, tmpStart, tmpEnd, tmpCN])
    return new_cnvList


def resolveConflictCNVs(cnvList):
    """Find and remove conflicted CNVs (overlapping DEL+DUP in same region)"""
    new_cnvList = []
    for x in cnvList:
        cnvLen = int(x[2]) - int(x[1])
        flag = 0
        for y in cnvList:
            if x != y:
                if testOverlapped(x[:3], y[:3]):
                    flag = 1
        if flag == 0 and cnvLen >= 100:
            new_cnvList.append(x)
    return new_cnvList


PRECISE_TOOLS = {'smoove', 'delly'}


def mergeCNVFromTools(cnvList, min_threshold=0.75, max_threshold=0.95, length_ratio_limit=None):
    """ Merge CNV results from different tools
    >>> mergeCNVFromTools([['chr1',0,1000,0,'smoove'], ['chr1',100,1000,0,'delly'], ['chr2',0,1000,0,'cnvkit']])
    [['chr1', 0, 1000, 0, 'delly,smoove', 2, 90.0], ['chr2', 0, 1000, 0, 'cnvkit', 1, 0.0]]
    """
    cnvs2 = cnvList[:]
    mergedCnvs = []
    while cnvList:
        cnv1 = cnvList.pop(0)
        cnvs2.pop(0)
        cnvs3 = cnvs2[:]
        count = 0
        accumLen = 0
        tmpCN = [cnv1[3]]
        for i, cnv2 in enumerate(cnvs3):
            overlap, prop1, prop2 = calculateOverlapProp(cnv1[:4], cnv2[:4])
            if min(prop1, prop2) > min_threshold or max(prop1, prop2) > max_threshold:
                # Length ratio check for cross-sample merging
                if length_ratio_limit is not None:
                    len1 = cnv1[2] - cnv1[1]
                    len2 = cnv2[2] - cnv2[1]
                    if min(len1, len2) / max(len1, len2) < length_ratio_limit:
                        continue
                # Only expand boundary if no precise tool (smoove/delly) has set it,
                # or the new tool is also precise
                current_tools = set(cnv1[-1].split(','))
                has_precise = bool(PRECISE_TOOLS & current_tools)
                new_is_precise = cnv2[-1] in PRECISE_TOOLS
                if not has_precise or new_is_precise:
                    cnv1[1:3] = [min(cnv1[1], cnv2[1]), max(cnv1[2], cnv2[2])]
                cnv1[-1] = ','.join([cnv1[-1], cnv2[-1]])
                cnvs2.pop(i - count)
                accumLen += overlap
                count += 1
                tmpCN.append(cnv2[3])

        tools = set(cnv1[-1].split(','))
        cnv1[-1] = ",".join(tools)
        cnv1.append(len(tools))
        cnv1[3] = max(tmpCN, key=tmpCN.count)
        accumFold = round(accumLen * 100 / (int(cnv1[2]) - int(cnv1[1])), 1)
        cnv1.append(accumFold)
        mergedCnvs.append(cnv1)
        cnvList = cnvs2[:]
    return mergedCnvs


if __name__ == '__main__':
    import doctest
    doctest.testmod()
