import argparse
import csv
import numpy as np

from ORForise.src.ORForise.utils import *  # local file

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome_to_compare', default='', help='Which genome to analyse?')
args = parser.parse_args()


def start_Codon_Count(start_Codons):
    atg, gtg, ttg, att, ctg, other = 0, 0, 0, 0, 0, 0
    other_Starts = []
    for start in start_Codons:
        if start == 'ATG':
            atg += 1
        elif start == 'GTG':
            gtg += 1
        elif start == 'TTG':
            ttg += 1
        elif start == 'ATT':
            att += 1
        elif start == 'CTG':
            ctg += 1
        else:
            other += 1
            other_Starts.append(start)
    atg_P = format(100 * atg / len(start_Codons), '.2f')
    gtg_P = format(100 * gtg / len(start_Codons), '.2f')
    ttg_P = format(100 * ttg / len(start_Codons), '.2f')
    att_P = format(100 * att / len(start_Codons), '.2f')
    ctg_P = format(100 * ctg / len(start_Codons), '.2f')
    other_Start_P = format(100 * other / len(start_Codons), '.2f')
    return atg_P, gtg_P, ttg_P, att_P, ctg_P, other_Start_P, other_Starts


def stop_Codon_Count(stop_Codons):
    tag, taa, tga, other = 0, 0, 0, 0
    other_Stops = []
    for stop in stop_Codons:
        stop
        if stop == 'TAG':
            tag += 1
        elif stop == 'TAA':
            taa += 1
        elif stop == 'TGA':
            tga += 1
        else:
            other += 1
            other_Stops.append(stop)
    tag_p = format(100 * tag / len(stop_Codons), '.2f')
    taa_p = format(100 * taa / len(stop_Codons), '.2f')
    tga_p = format(100 * tga / len(stop_Codons), '.2f')
    other_Stop_P = format(100 * other / len(stop_Codons), '.2f')
    return tag_p, taa_p, tga_p, other_Stop_P, other_Stops


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
    gc_content = (g + c) * 100 / (a + t + g + c + n)
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


def genome_Metrics(genome_to_compare):
    genome_Seq = ""
    with open('../Genomes/' + genome_to_compare + '.fa', 'r') as genome:
        for line in genome:
            line = line.replace("\n", "")
            if not line.startswith('>'):
                genome_Seq += str(line)

    genome_N_Per, genome_GC = gc_count(genome_Seq)

    genome_Rev = revCompIterative(genome_Seq)
    genome_Size = len(genome_Seq)
    coding_Regions = np.zeros((genome_Size), dtype=np.int)
    non_Coding_Regions = np.zeros((genome_Size), dtype=np.int)
    all_gene_Regions = np.zeros((genome_Size), dtype=np.int)
    protein_coding_genes = collections.OrderedDict()
    non_protein_coding_genes = collections.OrderedDict()
    strands = collections.defaultdict(int)
    lengths_PCG, gene_Pos_Olap, gene_Neg_Olap, short_PCGs, pcg_GC = [], [], [], [], []
    prev_Gene_Stop, count, nc_Count, pos_Strand, neg_Strand = 0, 0, 0, 0, 0
    prev_Gene_Overlapped = False
    with open('../Genomes/' + genome_to_compare + '.gff', 'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "CDS" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    length = stop - start
                    all_gene_Regions[start - 1:stop] = [1]
                    strand = line[6]
                    strands[strand] += 1
                    lengths_PCG.append(length)
                    coding_Regions[start - 1:stop] = [1]
                    gene = str(start) + ',' + str(stop) + ',' + strand
                    protein_coding_genes.update({count: gene})
                    if '+' in strand:
                        seq = genome_Seq[start - 1:stop]
                        pos_Strand += 1
                    elif '-' in strand:
                        r_Start = genome_Size - stop
                        r_Stop = genome_Size - start
                        seq = genome_Rev[r_Start:r_Stop + 1]
                        neg_Strand += 1
                    if length < SHORT_ORF_LENGTH:
                        short_PCGs.append(gene)
                    n_per, gc = gc_count(seq)
                    pcg_GC.append(gc)
                    ### Calculate overlapping ORFs -
                    if prev_Gene_Stop > start:
                        if '+' in strand:
                            gene_Pos_Olap.append(prev_Gene_Stop - start)
                        elif '-' in strand:
                            gene_Neg_Olap.append(prev_Gene_Stop - start)
                        prev_Gene_Overlapped = True
                    elif prev_Gene_Stop < start:
                        if prev_Gene_Overlapped == True:
                            if '+' in strand:
                                gene_Pos_Olap.append(0)
                            elif '-' in strand:
                                gene_Neg_Olap.append(0)
                        prev_Gene_Overlapped = False
                    prev_Gene_Stop = stop
                    count += 1
                elif "ID=gene" in line[8]:
                    gene_Info = line[8]
                    if "biotype=protein_coding" not in gene_Info:
                        start = int(line[3])
                        stop = int(line[4])
                        strand = line[6]
                        gene = str(start) + ',' + str(stop) + ',' + strand
                        all_gene_Regions[start - 1:stop] = [1]
                        non_Coding_Regions[start - 1:stop] = [1]
                        non_protein_coding_genes.update({nc_Count: gene})
                        nc_Count += 1

            except IndexError:
                continue

        if prev_Gene_Overlapped == True:  # If last has a prev overlap, count it
            if '+' in strand:
                gene_Pos_Olap.append(0)
            elif '-' in strand:
                gene_Neg_Olap.append(0)

    median_PCG = np.median(lengths_PCG)
    gene_Overlaps = gene_Neg_Olap + gene_Pos_Olap
    median_PCG_Olap = np.median(gene_Overlaps)
    longest_Olap = max(gene_Overlaps)
    coding_Percentage = 100 * float(np.count_nonzero(coding_Regions)) / float(genome_Size)
    non_coding_Percentage = 100 * float(np.count_nonzero(non_Coding_Regions)) / float(genome_Size)
    all_gene_Percentage = 100 * float(np.count_nonzero(all_gene_Regions)) / float(genome_Size)
    start_Codons, stop_Codons = [], []
    for gene in protein_coding_genes.values():
        start = int(gene.split(',')[0])
        stop = int(gene.split(',')[1])
        strand = gene.split(',')[2]

        if '-' in strand:
            r_start = genome_Size - stop
            r_stop = genome_Size - start

            start_Codons.append(genome_Rev[r_start:r_start + 3])
            stop_Codons.append(genome_Rev[r_stop - 2:r_stop + 1])

        elif '+' in strand:
            start_Codons.append(genome_Seq[start - 1:start - 1 + 3])
            stop_Codons.append(genome_Seq[stop - 3:stop - 1 + 1])

    atg_P, gtg_P, ttg_P, att_P, ctg_P, other_Start_P, other_Starts = start_Codon_Count(start_Codons)
    tag_P, taa_P, tga_P, other_Stop_P, other_Stops = stop_Codon_Count(stop_Codons)

    output = ("Number of Protein Coding Genes in " + genome_to_compare + " : " + str(
        len(lengths_PCG)) + " ,Median Length of PCGs: " + str(median_PCG) + ", Min Length of PCGs: " + str(
        min(lengths_PCG)) + ", Max Length of PCGs: " + str(max(lengths_PCG)) +
              ", Number of PCGs on Pos Strand: " + str(strands['+']) + ", Number of PCGs on Neg Strand: " + str(
                strands['-']) + "\nGenome-Wide GC Content: " + format(np.median(genome_GC),
                                                                      '.2f') + "Median GC of PCGs: " + format(
                np.median(pcg_GC), '.2f') + ", Number of Overlapping PCGs: " + str(len(gene_Overlaps)) +
              ", Longest PCG Overlap: " + str(longest_Olap) + ", Median PCG Overlap: " + str(
                median_PCG_Olap) + ", Number of PCGs less than 100 amino acids: " + str(len(short_PCGs)) +

              '\nPercentage of Genome which is Protein Coding: ' + format(coding_Percentage,
                                                                          '.2f') + ', Number of Non-PCGs: ' + str(
                len(non_protein_coding_genes)) + ', Percentage of Genome Non-PCG: ' + format(non_coding_Percentage,
                                                                                             '.2f') +
              ', Percentage of All Genes in Genome: ' + format(all_gene_Percentage, '.2f') +

              '\nPercentage of Genes starting with ATG: ' + atg_P +
              '\nPercentage of Genes starting with GTG: ' + gtg_P +
              '\nPercentage of Genes starting with TTG: ' + ttg_P +
              '\nPercentage of Genes starting with ATT: ' + att_P +
              '\nPercentage of Genes starting with CTG: ' + ctg_P +
              '\nPercentage of Genes starting with Alternative Start Codon: ' + other_Start_P +

              '\nPercentage of Genes ending with TAG: ' + tag_P +
              '\nPercentage of Genes ending with TAA: ' + taa_P +
              '\nPercentage of Genes ending with TGA: ' + tga_P +
              '\nPercentage of Genes ending with Alternative Stop Codon: ' + other_Stop_P)

    with open('../Genomes/' + genome_to_compare + '_metrics.csv', 'w') as out_file:
        out = csv.writer(out_file, delimiter=',')
        out.writerow(['Genome Metrics:'])
        out.writerow([output])

    print(output)


if __name__ == "__main__":
    genome_Metrics(**vars(args))
