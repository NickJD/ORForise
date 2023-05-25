#!/usr/bin/env python3
import collections

# Constants
SHORT_ORF_LENGTH = 300
MIN_COVERAGE = 75
ORForise_Version = 'v1.4.1'


def revCompIterative(watson):  # Gets Reverse Complement
    return watson.upper()[::-1].translate(str.maketrans("ATCGRYKMVBHD","TAGCYRMKBVDH"))


def sortORFs(tool_ORFs):  # Will only sort by given start position
    tool_ORFs_Sorted = sorted(tool_ORFs.items(), key=lambda v: int(v[0].split(",")[0]))
    tool_ORFs_Sorted = collections.OrderedDict(tool_ORFs_Sorted)
    return tool_ORFs_Sorted


def sortGenes(Genes):  # Will sort by given start position and then rearrange for given stop
    Genes_Sorted_list = sorted(Genes.values(), key=lambda v: int(v[0]))
    Genes_Sorted = []
    for idx,gene in enumerate(Genes_Sorted_list):
        Genes_Sorted.append([idx,gene])
    Genes_Sorted = collections.OrderedDict(Genes_Sorted)
    prev_stop = 0
    for pos, detail in Genes_Sorted.items():
        if detail[1] < prev_stop:
            Genes_Sorted[pos], Genes_Sorted[pos-1] = Genes_Sorted[pos-1], Genes_Sorted[pos]
        prev_stop = detail[1]
    return Genes_Sorted