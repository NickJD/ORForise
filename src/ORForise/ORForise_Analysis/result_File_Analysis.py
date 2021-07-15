import copy

import argparse
import collections

from ORForise.src.ORForise.utils import *

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tool', required=True, help='Which tool to compare?')
parser.add_argument('-p', '--parameters', required=False, help='Optional parameters for prediction tool.')
parser.add_argument('-g', '--genome', required=True, help='Which genome to analyse?')
args = parser.parse_args()


def gc_count(dna):
    c = 0
    a = 0
    g = 0
    t = 0
    n = 0
    for i in dna:
        if "C" in i:
            c += 1
        elif "G" in i:
            g += 1
        elif "A" in i:
            a += 1
        elif "T" in i:
            t += 1
        elif "N" in i:
            n += 1
    gc_content = format((g + c) * 100 / (a + t + g + c + n), '.2f')
    n_per = n * 100 / (a + t + g + c + n)
    return n_per, gc_content


def start_Codon_Count(orfs):
    atg, gtg, ttg, att, ctg, other = 0, 0, 0, 0, 0, 0
    other_Starts = []
    for orf in orfs.values():
        codon = orf[-2]
        if codon == 'ATG':
            atg += 1
        elif codon == 'GTG':
            gtg += 1
        elif codon == 'TTG':
            ttg += 1
        elif codon == 'ATT':
            att += 1
        elif codon == 'CTG':
            ctg += 1
        else:
            other += 1
            other_Starts.append(codon)
    atg_P = format(100 * atg / len(orfs), '.2f')
    gtg_P = format(100 * gtg / len(orfs), '.2f')
    ttg_P = format(100 * ttg / len(orfs), '.2f')
    att_P = format(100 * att / len(orfs), '.2f')
    ctg_P = format(100 * ctg / len(orfs), '.2f')
    other_Start_P = format(100 * other / len(orfs), '.2f')
    return atg_P, gtg_P, ttg_P, att_P, ctg_P, other_Start_P, other_Starts


def stop_Codon_Count(orfs):
    tag, taa, tga, other = 0, 0, 0, 0
    other_Stops = []
    for orf in orfs.values():
        codon = orf[-1]
        if codon == 'TAG':
            tag += 1
        elif codon == 'TAA':
            taa += 1
        elif codon == 'TGA':
            tga += 1
        else:
            other += 1
            other_Stops.append(codon)
    tag_p = format(100 * tag / len(orfs), '.2f')
    taa_p = format(100 * taa / len(orfs), '.2f')
    tga_p = format(100 * tga / len(orfs), '.2f')
    other_Stop_P = format(100 * other / len(orfs), '.2f')
    return tag_p, taa_p, tga_p, other_Stop_P, other_Stops


def detail_transfer(genes, missed_genes):
    for missed, m_details in missed_genes.items():
        try:
            details = genes[missed]
            gc = details[2]
            up_Overlap = details[3]
            down_Overlap = details[4]
            m_details.insert(2, gc)
            m_details.insert(3, up_Overlap)
            m_details.insert(4, down_Overlap)
        except KeyError:
            pass
    return missed_genes


def get_genome(genome):
    genome_Seq = ""
    with open('../Genomes/' + genome + '.fa', 'r') as genome:
        for line in genome:
            line = line.replace("\n", "")
            if not line.startswith('>'):
                genome_Seq += str(line)
    return genome_Seq


def missed_genes_in(genes_detected, missed_genes, results_in):
    # Missed Genes Read-In
    read = False
    for line in results_in:
        line = line.strip()
        if read == True:
            if line.startswith('>'):
                entry = line.split('_')
                entry = entry[1] + '_' + entry[2]
                strand = entry[-1]
                if int(strand) <= 2:
                    strand = '+'
                else:
                    strand = '-'
            elif len(line.strip()) > 0:
                startCodon = line[0:3]
                stopCodon = line[-3:]
                length = len(line)
                missed_genes.update({entry: [line, strand, length, startCodon, stopCodon]})

        if line.startswith('Undetected_Genes:'):
            read = True
        if read == True and not line:
            break
    list_Missed = list(missed_genes.keys())
    for key in list_Missed:
        # ### printed out to confirm figure lengths
        # start = key.split('_')[0]
        # stop = key.split('_')[1]
        if key in genes_detected:
            del genes_detected[key]

    return missed_genes, genes_detected


def partial_matches_in(partial_matches, results_in):
    # partial Genes Read-In
    read = False
    prev = ''
    for line in results_in:
        line = line.strip()
        if read == True:
            if line.startswith('Gene:'):
                line = line.replace('Gene:', '')
                entry = line.split('_')
                g_Pos = entry[0] + '_' + entry[1]
                strand = entry[2]
                prev = 'Gene'
            elif line.startswith('ORF:'):
                line = line.replace('ORF:', '')
                entry = line.split('_')
                o_Pos = entry[0] + '_' + entry[1]
                prev = 'ORF'
            elif line:
                if 'Gene' in prev:
                    g_Seq = line.strip()
                    g_length = len(g_Seq)
                elif 'ORF' in prev:
                    o_Seq = line.strip()
                    orf_length = len(o_Seq)
            elif not line:
                partial_matches.update({g_Pos: [strand, g_length, g_Seq, o_Pos, orf_length, o_Seq]})
        if line.startswith('Partial_Gene_Hits:'):
            read = True

    return partial_matches


def unmatched_ORFs_in(unmatched_ORFs, results_file):
    # Unmatched ORFs Read-In
    read = False
    for line in results_file:
        line = line.strip()
        if read == True:
            if line.startswith('>'):
                line = line.replace('Gene:', '')
                entry = line.split('_')
                strand = entry[-1]
                o_Pos = entry[1] + '_' + entry[2]
                unmatched_ORFs.update({o_Pos: [strand, None, None, None, None]})
            elif line:
                o_Seq = line.strip()
                o_Length = len(o_Seq)
                startCodon = line[0:3]
                stopCodon = line[-3:]
                unmatched_ORFs.update({o_Pos: [strand, o_Length, o_Seq, startCodon, stopCodon]})
        if line.startswith('ORF_Without_Corresponding_Gene_in_Ensembl:'):
            read = True
        elif read == True and not line:
            unmatched_ORFs.update({o_Pos: [strand, o_Length, o_Seq, startCodon, stopCodon]})
            break

    return unmatched_ORFs


def genes_in(genome, genome_Seq, genome_Seq_Rev, genome_Size, genes):
    with open('../Genomes/' + genome + '.gff', 'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "CDS" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]
                    length = stop - start
                    gene = str(start) + '_' + str(stop)
                    if '+' in strand:
                        seq = genome_Seq[start - 1:stop]
                    elif '-' in strand:
                        r_Start = genome_Size - stop
                        r_Stop = genome_Size - start
                        seq = genome_Seq_Rev[r_Start:r_Stop + 1]
                    startCodon = seq[0:3]
                    stopCodon = seq[-3:]
                    genes.update({gene: [seq, strand, length, startCodon, stopCodon]})
            except IndexError:
                continue
    return genes


def extract_results(genome, results_file):
    genome_Seq = get_genome(genome)
    genome_Seq_Rev = revCompIterative(genome_Seq)
    genome_Size = len(genome_Seq)
    genes = collections.OrderedDict()
    partial_matches = collections.OrderedDict()
    missed_genes = collections.OrderedDict()
    unmatched_ORFs = collections.OrderedDict()

    genes = genes_in(genome, genome_Seq, genome_Seq_Rev, genome_Size, genes)
    genes_detected = copy.deepcopy(genes)
    missed_genes, genes_detected = missed_genes_in(genes_detected, missed_genes, results_file)
    results_file.seek(0, 0)  # Reset file position
    partial_matches = partial_matches_in(partial_matches, results_file)
    results_file.seek(0, 0)  # Reset file position
    unmatched_ORFs = unmatched_ORFs_in(unmatched_ORFs, results_file)
    results_file.seek(0, 0)  # Reset file position

    return genes, genes_detected, missed_genes, partial_matches, unmatched_ORFs


if __name__ == "__main__":
    options = parser.parse_args()
    parameters = options.parameters
    tool = options.tool
    genome = options.genome
    if parameters:
        results_file = open('../Tools/' + tool + '/' + tool + '_' + genome + '_' + parameters + '.csv')
    else:
        results_file = open('../Tools/' + tool + '/' + tool + '_' + genome + '.csv')

    genes, genes_detected, missed_genes, partial_matches, unmatched_ORFs = extract_results(genome, results_file)
    gene_Lengths, genes_detected_Lengths, missed_Lengths, partial_Lengths, unmatched_Lengths = [], [], [], [], []

    for pos in genes.keys():
        gene_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    for pos in genes_detected.keys():
        genes_detected_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    for pos in partial_matches.keys():
        partial_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    for pos in missed_genes.keys():
        missed_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    for pos in unmatched_ORFs.keys():
        unmatched_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))

    import numpy as np

    print(len(gene_Lengths))
    print(gene_Lengths)
    print(np.mean(gene_Lengths))
    print(len(partial_Lengths))
    print(partial_Lengths)
    print(len(missed_Lengths))
    print(missed_Lengths)
    print(len(unmatched_Lengths))
    print(unmatched_Lengths)
    print(np.mean(unmatched_Lengths))
