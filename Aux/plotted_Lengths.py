import argparse
import Tools.utils as utils
import collections
import numpy as np
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--tool', required=True, help='Which tool to plot?')
parser.add_argument('-g', '--genomes', required=True, help='Which genomes?')
args = parser.parse_args()


def get_Lengths(tool,genomes):
    genomes = genomes.split(',')
    all_Genomes = []
    for genome in genomes:

        gene_Pos = result_compare(genome)
        gene_Lengths = []
        missed_Pos = []
        missed_Lengths = []
        orf_Pos = []
        orf_Lengths = []
        partial_Pos = []
        partial_Lengths = []


        predictions_in = open('../Tools/' + tool + '/' + tool + '_' + genome + '.csv')
        missed = False
        orfs = False
        partials = False
        for line in predictions_in:
            line = line.strip()
            if line.startswith('Undetected_Genes:'):
                missed = True
            elif line.startswith('ORF_Without_Corresponding_Gene_in_Ensembl:'):
                missed = False
                orfs = True
            elif line.startswith('Partial_Gene_Hits:'):
                missed = False
                orfs = False
                partials = True

            if missed == True:
                if line.startswith('>'):
                    pos = line.split('_')[1] + '_' + line.split('_')[2]
                    print(pos)
                    try:
                        gene_Pos.remove(pos)
                    except ValueError:
                        print('Already Removed')
                    missed_Pos.append(pos)
            elif orfs == True:
                if line.startswith('>'):
                    pos = line.split('_')[1] + '_' + line.split('_')[2]
                    print(pos)
                    orf_Pos.append(pos)
            if partials == True:
                if line.startswith('Gene:'):
                    pos = line.split('_')[0].split(':')[1] + '_' + line.split('_')[1]
                    print(pos)
                    try:
                        gene_Pos.remove(pos)
                    except ValueError:
                        print('Already Removed')
                    partial_Pos.append(pos)


        for pos in gene_Pos:
            gene_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
        for pos in missed_Pos:
            missed_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
        for pos in orf_Pos:
            orf_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))
        for pos in partial_Pos:
            partial_Lengths.append(int(pos.split('_')[1]) - int(pos.split('_')[0]))


        genome = [gene_Lengths,missed_Lengths,orf_Lengths,partial_Lengths]
        all_Genomes.append(genome)

    return all_Genomes







def result_compare(genome):
    gene_Pos = []
    with open('../Genomes/' + genome + '.gff', 'r') as genome_file:
        for line in genome_file:
            if line.startswith('Chromosome	ena	CDS'):
                line = line.split()
                pos = line[3]+'_'+line[4]
                gene_Pos.append(pos)
    return  gene_Pos




if __name__ == "__main__":


    all_Genomes = get_Lengths(**vars(args))

    for genome in all_Genomes:
        for set in genome:
            print(set)





    ## plot

    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    bp0 = plt.boxplot(all_Genomes[0], positions=[0.1,0.4,0.7,1], sym='', widths=0.4,vert=False)
    bp1 = plt.boxplot(all_Genomes[1], positions=[1.8,2.1,2.4,2.7], sym='', widths=0.4,vert=False)
    bp2 = plt.boxplot(all_Genomes[2], positions=[3.8,4.1,4.4,4.7],sym='', widths=0.4, vert=False)
    bp3 = plt.boxplot(all_Genomes[3], positions=[5.8,6.1,6.4,6.7],sym='', widths=0.4, vert=False)
    bp4 = plt.boxplot(all_Genomes[4], positions=[7.8,8.1,8.4,8.7],sym='', widths=0.4, vert=False)
    bp5 = plt.boxplot(all_Genomes[5], positions=[9.8,10.1,10.4,10.7],sym='', widths=0.4, vert=False)

    ticks = ['E-coli', 'Staph', 'Bacillus', 'Myco', 'Caul', 'Pseudo']
    plt.plot([], c='#D7191C', label='Missed')
    plt.plot([], c='#2C7BB6', label='Covered')
    plt.yticks(range(0, len(ticks) * 2, 2), ticks)

    plt.tight_layout()
    plt.show()




