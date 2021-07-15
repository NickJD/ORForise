import argparse
import collections
import numpy as np

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


def revCompIterative(watson):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev:
        crick += complements[nt]
    return crick


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


def missed_gene_read(results_file, missed_genes):
    # Missed Genes Read-In
    read = False
    for line in results_file:
        if line.startswith('ORFs_Without_Corresponding_Gene_In_Ensembl_Metrics:'):
            break
        line = line.strip()
        if read == True:
            if line.startswith('>'):
                entry = line.split('_')
                entry = entry[1] + '_' + entry[2]
            elif len(line.strip()) > 0:
                startCodon = line[0:3]
                stopCodon = line[-3:]
                length = len(line)
                missed_genes.update({entry: [line, length, startCodon, stopCodon]})
        if line.startswith('Undetected_Genes:'):
            read = True

    return missed_genes


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


def result_compare(genome, results_file):
    genome_Seq = ""
    with open('../Genomes/' + genome + '.fa', 'r') as genome_file:
        for line in genome_file:
            line = line.replace("\n", "")
            if ">" not in line:
                genome_Seq += str(line)

    missed_genes = collections.OrderedDict()
    missed_genes = missed_gene_read(results_file, missed_genes)
    list_MG = list(missed_genes.keys())
    # Analysis
    genome_Rev = revCompIterative(genome_Seq)
    genome_Size = len(genome_Seq)
    genes = collections.OrderedDict()
    count = 0
    prev_Stop = 0
    ### Record Missed and Detected Gene metrics
    genes_strand = collections.defaultdict(int)
    genes_Missed_strand = collections.defaultdict(int)
    short_PCGs, short_Missed_PCGs, pcg_GC, pcg_Missed_GC, lengths_PCG, lengths_Missed_PCG, genes_Overlap, genes_Missed_Overlap = [], [], [], [], [], [], [], []
    with open('../Genomes/' + genome + '.gff', 'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "CDS" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]
                    gene = str(start) + ',' + str(stop) + ',' + strand
                    if strand == '-':
                        r_Start = genome_Size - stop
                        r_Stop = genome_Size - start
                        seq = (genome_Rev[r_Start:r_Stop + 1])
                    elif strand == '+':
                        seq = (genome_Seq[start - 1:stop])
                    startCodon = seq[0:3]
                    stopCodon = seq[-3:]
                    length = stop - start
                    n_per, gc = gc_count(seq)
                    pos = str(start) + '_' + str(stop)
                    if pos in list_MG:
                        genes_Missed_strand[strand] += 1
                        pcg_Missed_GC.append(float(gc))
                        lengths_Missed_PCG.append(length)
                        if length < SHORT_ORF_LENGTH:
                            short_Missed_PCGs.append(gene)
                        if prev_Stop > start:
                            overlap = prev_Stop - start
                            genes_Missed_Overlap.append(overlap)
                        else:
                            overlap = 0
                    elif pos not in list_MG:
                        genes_strand[strand] += 1
                        pcg_GC.append(float(gc))
                        lengths_PCG.append(length)
                        if length < SHORT_ORF_LENGTH:
                            short_PCGs.append(gene)
                        if prev_Stop > start:
                            overlap = prev_Stop - start
                            genes_Overlap.append(overlap)
                        else:
                            overlap = 0

                    count += 1
                    prev_Stop = stop
                    pos = str(start) + '_' + str(stop)
                    if genes:
                        prev_details = genes[prev_pos]
                        prev_details.insert(4, overlap)
                        genes.update({prev_pos: prev_details})

                    genes.update({pos: [strand, length, gc, overlap, seq, startCodon, stopCodon]})
                    prev_pos = pos

            except IndexError:
                continue

    missed_genes = detail_transfer(genes, missed_genes)

    missed_lengths = []
    gene_lengths = []
    for key in list_MG:
        ### printed out to confirm figure lengths
        start = key.split('_')[0]
        stop = key.split('_')[1]
        m_length = int(stop) - int(start)
        missed_lengths.append(m_length)
        if key in genes:
            del genes[key]

    ## Printed out for Figure of gene lengths
    for key in genes.keys():
        start = key.split('_')[0]
        stop = key.split('_')[1]
        g_length = int(stop) - int(start)
        gene_lengths.append(g_length)

    print("Number of Genes Missed:" + str(len(missed_lengths)) + '\nLengths of Genes Missed:\n' + str(missed_lengths))

    median_PCG = np.median(lengths_PCG)
    median_PCG_Olap = np.median(genes_Overlap)
    longest_Olap = max(genes_Overlap)
    num_overlaps = len(genes_Overlap)
    gc_median = format(np.median(pcg_GC), '.2f')
    num_Short_PCGs = len(short_PCGs)

    atg_P, gtg_P, ttg_P, att_P, ctg_P, other_Starts_P, other_Starts = start_Codon_Count(genes)
    tag_P, taa_P, tga_P, other_Stops_P, other_Stops = stop_Codon_Count(genes)
    m_atg_P, m_gtg_P, m_ttg_P, m_att_P, m_ctg_P, m_other_Starts_P, m_other_Starts = start_Codon_Count(missed_genes)
    m_tag_P, m_taa_P, m_tga_P, m_other_Stops_P, m_other_Stops = stop_Codon_Count(missed_genes)

    output = ("Number of Missed Protein Coding Genes in " + str(genome) + " : " + str(
        len(lengths_PCG)) + ", Median Length of PCGs: " +
              str(median_PCG) + ", Min Length of PCGs: " + str(min(lengths_PCG)) + ", Max Length of PCGs: " + str(
                max(lengths_PCG)) +
              ", Number of PCGs on Pos Strand: " + str(
                genes_Missed_strand['+']) + ", Number of PCGs on Neg Strand: " + str(genes_Missed_strand['-']) +
              ", Median GC of PCGs: " + str(gc_median) + ", Number of Overlapping PCGs: " + str(num_overlaps) +
              ", Longest PCG Overlap: " + str(longest_Olap) + ", Median PCG Overlap: " + str(median_PCG_Olap) +
              ", Number of PCGs less than 100nt: " + str(num_Short_PCGs) +

              '\nPercentage of Genes starting with ATG - Annotation/Missed: ' + atg_P + ' ' + m_atg_P +
              '\nPercentage of Genes starting with GTG - Annotation/Missed: ' + gtg_P + ' ' + m_gtg_P +
              '\nPercentage of Genes starting with TTG - Annotation/Missed: ' + ttg_P + ' ' + m_ttg_P +
              '\nPercentage of Genes starting with ATT - Annotation/Missed: ' + att_P + ' ' + m_att_P +
              '\nPercentage of Genes starting with CTG - Annotation/Missed: ' + ctg_P + ' ' + m_ctg_P +
              '\nPercentage of Genes starting with Alternative Start Codon - Annotation/Missed: ' + other_Starts_P + ' ' + m_other_Stops_P +
              '\nPercentage of Genes ending with TAG - Annotation/Missed: ' + tag_P + ' ' + m_tag_P +
              '\nPercentage of Genes ending with TAA - Annotation/Missed: ' + taa_P + ' ' + m_taa_P +
              '\nPercentage of Genes ending with TGA - Annotation/Missed: ' + tga_P + ' ' + m_tga_P +
              '\nPercentage of Genes ending with Alternative Stop Codon - Annotation/Missed: ' + other_Stops_P + ' ' + m_other_Stops_P)

    print(output)


if __name__ == "__main__":
    options = parser.parse_args()
    parameters = options.parameters
    tool = options.tool
    genome = options.genome
    if parameters:
        results_file = open('../Tools/' + tool + '/' + tool + '_' + genome + '_' + parameters + '.csv')
    else:
        results_file = open('../Tools/' + tool + '/' + tool + '_' + genome + '.csv')
    result_compare(genome, results_file)
