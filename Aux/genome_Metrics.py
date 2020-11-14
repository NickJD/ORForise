import argparse
import collections
import csv

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome_to_compare', default='', help='Which genome to analyse?')
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
            t+=1
        elif "N" in i:
            n+=1
    gc_content = (g + c) * 100 / (a + t + g + c + n)
    n_per = n * 100 / (a + t + g + c + n)
    return n_per,gc_content

def revCompIterative(watson):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev:
        crick += complements[nt]
    return crick

def genome_Metrics(genome_to_compare):
    genome = ""
    with open('../genomes/' + genome_to_compare + '.fa', 'r') as genome_file:
        for line in genome_file:
            line = line.replace("\n", "")
            if ">" not in line:
                genome += str(line)

    genome_Rev = revCompIterative(genome)
    genome_Size = len(genome)
    coding_Regions = np.zeros((genome_Size), dtype=np.int)
    non_Coding_Regions = np.zeros((genome_Size), dtype=np.int)
    all_gene_Regions = np.zeros((genome_Size), dtype=np.int)
    protein_coding_genes = collections.OrderedDict()
    non_protein_coding_genes = []
    lengths_PCG = []
    gene_Overlaps = []
    count = 0
    prev_Stop = 0
    strands = collections.defaultdict(int)
    short_PCGs = []
    pcg_GC = []
    with open('../genomes/' + genome_to_compare + '.gff', 'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "CDS" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    all_gene_Regions[start:stop] = [1]
                    strand = line[6]
                    strands[strand] +=1
                    gene = str(start) + ',' + str(stop) + ',' + strand
                    if '+' in strand:
                        seq = genome[start-1:stop]
                    elif '-' in strand:
                        seq = genome_Rev[start:stop]
                    length = stop-start
                    if length < 100:
                        short_PCGs.append(gene)
                        print(line)
                    n_per, gc = gc_count(seq)
                    pcg_GC.append(gc)
                    lengths_PCG.append(length)
                    coding_Regions[start-1:stop] = [1]
                    protein_coding_genes.update({count: gene})
                    if prev_Stop > start:
                        overlap = prev_Stop - start
                        gene_Overlaps.append(overlap)
                    count += 1
                    prev_Stop = stop

                elif "ID=gene" in line[8]:
                    gene_Info = line[8]
                    if "biotype=protein_coding" not in gene_Info:
                        start = int(line[3])
                        stop = int(line[4])
                        all_gene_Regions[start:stop] = [1]
                        non_Coding_Regions[start - 1:stop] = [1]
                        non_protein_coding_genes.append(gene)

            except IndexError:
                continue

    median_PCG = np.median(lengths_PCG)
    median_PCG_Olap = np.median(gene_Overlaps)
    longest_Olap = max(gene_Overlaps)
    coding_Percentage = 100 * float(np.count_nonzero(coding_Regions))/float(genome_Size)
    non_coding_Percentage = 100 * float(np.count_nonzero(non_Coding_Regions))/float(genome_Size)
    all_gene_Percentage = 100 * float(np.count_nonzero(all_gene_Regions))/float(genome_Size)
    starts, stops = [],[]
    for gene in protein_coding_genes.values():
        start = int(gene.split(',')[0])
        stop = int(gene.split(',')[1])
        strand = gene.split(',')[2]

        if '-' in strand:
            r_start = genome_Size - stop
            r_stop = genome_Size - start

            starts.append(genome_Rev[r_start:r_start + 3])
            stops.append(genome_Rev[r_stop - 2:r_stop + 1])

        elif '+' in strand:
            starts.append(genome[start-1:start-1 + 3])
            stops.append(genome[stop - 3:stop - 1 + 1])

    ATG = 0
    GTG = 0
    TTG = 0
    ATT = 0
    CTG = 0
    Other = 0

    for codon in starts:

        if codon == 'ATG':
            ATG += 1
        elif codon == 'GTG':
            GTG += 1
        elif codon == 'TTG':
            TTG += 1
        elif codon == 'ATT':
            ATT += 1
        elif codon == 'CTG':
            CTG += 1
        else:
            Other += 1
            print(codon)

    ATG_Percentage = float((ATG) * float(100) / float(len(lengths_PCG)))
    GTG_Percentage = float((GTG) * float(100) / float(len(lengths_PCG)))
    TTG_Percentage = float((TTG) * float(100) / float(len(lengths_PCG)))
    ATT_Percentage = float((ATT) * float(100) / float(len(lengths_PCG)))
    CTG_Percentage = float((CTG) * float(100) / float(len(lengths_PCG)))
    other_Start_Percentage = float(Other) * float(100) / float(len(lengths_PCG))

    TAG = 0
    TAA = 0
    TGA = 0
    Other = 0
    count = 0
    for codon in stops:
        count +=1
        if codon == 'TAG':
            TAG += 1
        elif codon == 'TAA':
            TAA += 1
        elif codon == 'TGA':
            TGA += 1
        else:
            Other +=1

    TAG_Percentage = float((TAG) * float(100) / float(len(lengths_PCG)))
    TAA_Percentage = float((TAA) * float(100) / float(len(lengths_PCG)))
    TGA_Percentage = float((TGA) * float(100) / float(len(lengths_PCG)))
    other_Stop_Percentage = float(Other) * float(100) / float(len(lengths_PCG))

    output = ("Number of Protein Coding Genes in " + genome_to_compare + " : " + str(len(lengths_PCG)) + ", Median Length of PCGs: " + str(median_PCG) + ", Min Length of PCGs: " + str(min(lengths_PCG)) + ", Max Length of PCGs: " + str(max(lengths_PCG)) +
              ", Number of PCGs on Pos Strand: " + str(strands['+']) + ", Number of PCGs on Neg Strand: " + str(strands['-']) + ", Median GC of PCGs: " + format(np.median(pcg_GC), '.2f') + ", Number of Overlapping PCGs: " + str(len(gene_Overlaps)) +
               ", Longest PCG Overlap: " + str(longest_Olap) + ", Median PCG Overlap: " + str(median_PCG_Olap)  + ", Number of PCGs less than 100nt: " + str(len(short_PCGs)) +

            '\nPercentage of Genome which is Protein Coding: ' + format(coding_Percentage, '.2f') +', Number of Non-PCGs: ' + str(len(non_protein_coding_genes)) + ', Percentage of Genome Non-PCG: ' +format(non_coding_Percentage, '.2f') +
            ', Percentage of All Genes in Genome: ' + format(all_gene_Percentage,'.2f') +

            '\nPercentage of Genes starting with ATG: ' + format(ATG_Percentage,'.2f') +
            '\nPercentage of Genes starting with GTG: ' + format(GTG_Percentage, '.2f') +
            '\nPercentage of Genes starting with TTG: ' + format(TTG_Percentage,'.2f') +
            '\nPercentage of Genes starting with ATT: ' + format(ATT_Percentage, '.2f') +
            '\nPercentage of Genes starting with CTG: ' + format(CTG_Percentage,'.2f') +
            '\nPercentage of Genes starting with Alternative Start Codon: ' + format(other_Start_Percentage, '.2f') +

            '\nPercentage of Genes ending with TAG: ' + format(TAG_Percentage,'.2f') +
            '\nPercentage of Genes ending with TAA: ' + format(TAA_Percentage, '.2f') +
            '\nPercentage of Genes ending with TGA: ' + format(TGA_Percentage, '.2f') +
            '\nPercentage of Genes ending with Alternative Stop Codon: ' + format(other_Stop_Percentage, '.2f'))

    with open('../genomes/' + genome_to_compare + '_metrics.csv', 'w') as out_file:
        out = csv.writer(out_file, delimiter=',')
        out.writerow(['Genome Metrics:'])
        out.writerow([output])


if __name__ == "__main__":
    genome_Metrics(**vars(args))
