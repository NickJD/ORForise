import argparse
from PCG_Comparison.Tools.utils import  *
import collections
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--tool', required=True, help='Which tool to plot?')
parser.add_argument('-g', '--genome', required=True, help='Which genome?')
args = parser.parse_args()

def detail_transfer(genes,missed_genes):
    for missed,m_details in missed_genes.items():
        try:
            details = genes[missed]
            gc = details[2]
            up_Overlap = details[3]
            down_Overlap = details[4]
            m_details.insert(2,gc)
            m_details.insert(3,up_Overlap)
            m_details.insert(4,down_Overlap)
        except KeyError:
            pass
    return missed_genes


def get_genome(genome):
    genome_Seq = ""
    with open('../Genomes/'+genome+'.fa', 'r') as genome:
        for line in genome:
            line = line.replace("\n","")
            if not line.startswith('>'):
                genome_Seq += str(line)
    return genome_Seq


def missed_gene_read(tool,missed_genes,genome):
    #Missed Genes Read-In
    results_in = open('../Tools/'+tool+'/'+tool+'_'+genome+'.csv')
    read = False
    for line in results_in:
        line = line.strip()
        if read == True:
            if line.startswith('>'):
                entry = line.split('_')
                entry = entry[1]+'_'+entry[2]
            elif len(line.strip()) > 0:
                startCodon = line[0:3]
                stopCodon = line[-3:]
                length = len(line)
                missed_genes.update({entry:[line,length,startCodon,stopCodon]})

        if line.startswith('Undetected_Genes:'):
            read = True
        if read == True and not line:
            break

    return missed_genes



def result_compare(tool,genome):
    missed_genes = collections.OrderedDict()
    missed_genes = missed_gene_read(tool,missed_genes,genome)
    list_MG = list(missed_genes.keys())
    #Analysis
    genome_Seq = get_genome(genome)
    genome_Rev = revCompIterative(genome_Seq)
    genome_Size = len(genome_Seq)
    genes = collections.OrderedDict()
    count = 0
    prev_Stop = 0
    ### Record Missed and Detected Gene metrics
    genes_strand = collections.defaultdict(int)
    genes_Missed_strand = collections.defaultdict(int)
    short_PCGs,short_Missed_PCGs,pcg_GC,pcg_Missed_GC,lengths_PCG,lengths_Missed_PCG,genes_Overlap,genes_Missed_Overlap = [],[],[],[],[],[],[],[]
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
                    pos = str(start) + '_' + str(stop)

                    if pos in list_MG:
                        genes_Missed_strand[strand] += 1
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
                    pos = str(start)+'_'+str(stop)
                    # if genes:
                    #     prev_details = genes[prev_pos]
                    #     prev_details.insert(4,overlap)
                    #     genes.update({prev_pos:prev_details})

                    genes.update({pos:[strand,length,overlap,seq,startCodon,stopCodon]})
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


    return genes, missed_genes


def partial_gene_read(tool,genome,genes):
    #partial Genes Read-In
    partial_matches = collections.OrderedDict()
    orf_Lengths = []
    gene_Lengths = []
    results_in = open('../Tools/'+tool+'/'+tool+'_'+genome+'.csv')
    read = False
    prev = ''
    for line in results_in:
        line = line.strip()
        if read == True:
            if line.startswith('Gene:'):
                line = line.replace('Gene:','')
                entry = line.split('_')
                g_Pos = entry[0]+'_'+entry[1]
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
                    gene_Lengths.append(len(g_Seq))
                elif 'ORF' in prev:
                    o_Seq = line.strip()
                    orf_Lengths.append(len(o_Seq))
            elif not line:
                partial_matches.update({g_Pos:[strand,g_Seq,o_Pos,o_Seq]})
        if line.startswith('Partial_Gene_Hits:'):
            read = True

    list_PM = list(partial_matches.keys())
    for key in list_PM:
        ### printed out to confirm figure lengths
        start = key.split('_')[0]
        stop = key.split('_')[1]
        if key in genes:
            del genes[key]



    return genes,partial_matches


def get_Lengths(tool,genome):

    genes, missed_genes = result_compare(tool,genome)

    genes,partial_matches = partial_gene_read(tool,genome,genes)




    return genes, missed_genes, partial_matches


    # predictions_in = open('../Tools/' + tool + '/' + tool + '_' + genome + '.csv')
    # missed = False
    # orfs = False
    # partials = False
    # for line in predictions_in:
    #     line = line.strip()
    #     if line.startswith('Undetected_Genes:'):
    #         missed = True
    #     elif line.startswith('ORF_Without_Corresponding_Gene_in_Ensembl:'):
    #         missed = False
    #         orfs = True
    #     elif line.startswith('Partial_Gene_Hits:'):
    #         missed = False
    #         orfs = False
    #         partials = True
    #
    #     if missed == True:
    #         if line.startswith('>'):
    #             pos = line.split('_')[1] + '_' + line.split('_')[2]
    #             print(pos)
    #             try:
    #                 gene_Pos.remove(pos)
    #             except ValueError:
    #                 print('Already Removed')
    #             missed_Pos.append(pos)
    #     elif orfs == True:
    #         if line.startswith('>'):
    #             pos = line.split('_')[1] + '_' + line.split('_')[2]
    #             print(pos)
    #             orf_Pos.append(pos)
    #     if partials == True:
    #         if line.startswith('Gene:'):
    #             pos = line.split('_')[0].split(':')[1] + '_' + line.split('_')[1]
    #             print(pos)
    #             try:
    #                 gene_Pos.remove(pos)
    #             except ValueError:
    #                 print('Already Removed')
    #             partial_Pos.append(pos)
    #
    #
    # # for pos in gene_Pos:
    # #     gene_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    # for pos in missed_Pos:
    #     missed_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    # # for pos in orf_Pos:
    # #     orf_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    # for pos in partial_Pos:
    #     partial_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    #
    #
    # #genome = [gene_Lengths,missed_Lengths,orf_Lengths,partial_Lengths]
    # genome = [missed_Lengths,partial_Lengths]
    # all_Genomes.append(genome)









# def result_compare(genome):
#     gene_Pos = []
#     with open('../Genomes/' + genome + '.gff', 'r') as genome_file:
#         for line in genome_file:
#             if line.startswith('Chromosome	ena	CDS'):
#                 line = line.split()
#                 pos = line[3]+'_'+line[4]
#                 gene_Pos.append(pos)
#     return  gene_Pos




if __name__ == "__main__":


    genes,missed_genes,partial_matches = get_Lengths(**vars(args))

    gene_Lengths = []
    missed_Lengths = []
    partial_Lengths = []




    for pos in genes.keys():
        gene_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    for pos in partial_matches.keys():
        partial_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
    for pos in missed_genes.keys():
        missed_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))

    gene_set = set(genes.keys())
    partial_set = set(partial_matches.keys())
    missed_set = set(missed_genes.keys())

    intersection = set.intersection(gene_set,partial_set,missed_set)
    if intersection:
        sys.exit('Intersection Error:')

    gene_Starts = collections.defaultdict(int)

    for gene, data in genes.items():

        gene_start = data[4]
        gene_stop = data[5]
        gene_Starts[gene_start] +=1

    for start, number in gene_Starts.items():
        print(start+':'+str(number))

    print(len(gene_Lengths))
    print(gene_Lengths)
    print(len(partial_Lengths))
    print(partial_Lengths)
    print(len(missed_Lengths))
    print(missed_Lengths)


    ## plot

    # import seaborn as sns
    # import matplotlib.pyplot as plt
    # import numpy as np
    #
    # bp0 = plt.boxplot(all_Genomes[0], positions=[0.1,0.4,0.7,1], sym='', widths=0.4,vert=False)
    # bp1 = plt.boxplot(all_Genomes[1], positions=[1.8,2.1,2.4,2.7], sym='', widths=0.4,vert=False)
    # bp2 = plt.boxplot(all_Genomes[2], positions=[3.8,4.1,4.4,4.7],sym='', widths=0.4, vert=False)
    # bp3 = plt.boxplot(all_Genomes[3], positions=[5.8,6.1,6.4,6.7],sym='', widths=0.4, vert=False)
    # bp4 = plt.boxplot(all_Genomes[4], positions=[7.8,8.1,8.4,8.7],sym='', widths=0.4, vert=False)
    # bp5 = plt.boxplot(all_Genomes[5], positions=[9.8,10.1,10.4,10.7],sym='', widths=0.4, vert=False)
    #
    # ticks = ['E-coli', 'Staph', 'Bacillus', 'Myco', 'Caul', 'Pseudo']
    # plt.plot([], c='#D7191C', label='Missed')
    # plt.plot([], c='#2C7BB6', label='Covered')
    # plt.yticks(range(0, len(ticks) * 2, 2), ticks)
    #
    # plt.tight_layout()
    # plt.show()
    #
    #
    #

