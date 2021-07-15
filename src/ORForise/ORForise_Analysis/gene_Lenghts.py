import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome_to_compare', default='', help='Which genome to analyse?')
args = parser.parse_args()


def genome_Lengths(genome_to_compare):
    lengths = []
    with open('../Genomes/' + genome_to_compare + '.gff', 'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "CDS" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    length = stop - start
                    lengths.append(length)
            except IndexError:
                # print(line)
                continue
    print("Number of Genes: " + str(len(lengths)) +
          '\tMedian Length of Genes: ' + str(np.median(lengths)) + '\nGenes Lengths:\n' + str(lengths))


if __name__ == "__main__":
    genome_Lengths(**vars(args))
