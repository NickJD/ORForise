import argparse
import collections
import csv
from PCG_Comparison.Tools.utils import * # local file
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome_to_compare', default='', help='Which genome to analyse?')
args = parser.parse_args()


def genome_Metrics(genome_to_compare):
    genome_Seq = ""
    with open('../Genomes/'+genome_to_compare+'.fa', 'r') as genome:
        for line in genome:
            line = line.replace("\n","")
            if not line.startswith('>'):
                genome_Seq += str(line)

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
                print(line)
                continue


    print(lengths)
    print(len(lengths))



if __name__ == "__main__":
    genome_Metrics(**vars(args))
