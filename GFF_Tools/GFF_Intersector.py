import argparse
import collections
from importlib import import_module
import numpy as np

import sys
sys.path.append('/home/nick/Git/')
sys.path.append('.')
import os
os.chdir('..')

from ORForise.utils import sortORFs

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome', required=True, help='Which genome to use as Gold Standard?')
parser.add_argument('-t', '--tool', required=True, help='While tool to add to Gold Standard Annotation?')
parser.add_argument('-p', '--parameters', required=False, help='Optional parameters for prediction tool.')
parser.add_argument('-cov', '--coverage', required=False, action='store', dest='coverage', default='100', type=int,
                        help='ORF coverage of Gene in % ')
parser.add_argument('-o', '--out_file', required=True, help='output filename')
args = parser.parse_args()


def gff_writer(genome,tool,outfile,new_ORFs):
    write_out = open(outfile, 'w')
    for pos in new_ORFs:
        pos_ = pos.split(',')
        start = pos_[0]
        stop = pos_[-1]
        type = 'original'
        entry = (genome + '\t' + type+ '\tORF\t' + start + '\t' + stop + '\t.\t' + '\n')
        write_out.write(entry)

def comparator(genome,tool,parameters,coverage,out_file): # Only works for single contig genome
    genome_Seq = ""
    with open('Genomes/'+genome+'.fa', 'r') as genome_fasta:
        for line in genome_fasta:
            line = line.replace("\n","")
            if not line.startswith('>'):
                genome_Seq += str(line)
    genome_Size = len(genome_Seq)
    gene_Nuc_Array = np.zeros((genome_Size), dtype=np.int)
    ###########################################
    gold_standard = collections.OrderedDict() # Order is important
    with open('Genomes/'+genome+'.gff','r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "CDS" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]
                    gene_Nuc_Array[start - 1:stop] = [1]  # Changing all between the two positions to 1's
                    pos = str(start) + ',' + str(stop)
                    gold_standard.update({pos:strand})
            except IndexError:
                pass

    tool_predictions = import_module('Tools.'+tool+'.'+tool)
    tool_predictions = getattr(tool_predictions,tool)
    orfs = tool_predictions(genome,parameters,genome_Seq)
    genes_To_Keep = []


    if coverage == 100:
        for orf in orfs.keys():
            o_Start = int(orf.split(',')[0])
            o_Stop = int(orf.split(',')[1])
            try:
                if gold_standard[str(o_Start)+','+str(o_Stop)]:
                    genes_To_Keep.append(str(o_Start)+','+str(o_Stop)) #  o_ and g_ would be the same here
            except KeyError:
                continue
    else:
        for orf, data in orfs.items():
            o_Start = int(orf.split(',')[0])
            o_Stop = int(orf.split(',')[1])
            o_Strand = data[0]
            orf_Set = set(range(int(o_Start), int(o_Stop) + 1))
            for gene, g_data in gold_standard.items():  # Very ineffecient
                g_Start = int(gene.split(',')[0])
                g_Stop = int(gene.split(',')[1])
                g_Strand = g_data[0]
                gene_Set = set(range(int(g_Start), int(g_Stop) + 1))
                overlap = len(orf_Set.intersection(gene_Set))
                cov = 100 * float(overlap) / float(len(gene_Set))
                if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and cov >= coverage:
                    genes_To_Keep.append(str(g_Start) + ',' + str(g_Stop))
                if g_Start > o_Stop:
                    break
    #########################################################
    gff_writer(genome,tool,out_file,genes_To_Keep)

if __name__ == "__main__":
    comparator(**vars(args))

    print("Complete")
