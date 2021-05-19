import argparse
import collections
from importlib import import_module
import numpy as np
from datetime import date
from utils import sortORFs

parser = argparse.ArgumentParser()
parser.add_argument('-g', action='store', dest='genome', required=True,
                    help='Which genome to use as Gold Standard?')
parser.add_argument('-t', action='store', dest='tool', required=True,
                    help='While tool to add to Gold Standard Annotation?')
parser.add_argument('-p', action='store', dest='parameters', required=False,
                    help='Optional parameters for prediction tool.')
parser.add_argument('-olap', action='store', dest='overlap', default=50, type=int, required=False,
                    help='maximum overlap between Gene and ORF')
parser.add_argument('-o', action='store', dest='output_file',  required=True,
                    help='output filename')
args = parser.parse_args()


def gff_writer(genome,tool,output_file,new_ORFs,genome_gff):
    write_out = open(output_file, 'w')
    write_out.write("##gff-version\t3\n#\tGFF Adder\n#\tRun Date:" + str(date.today()) + '\n')
    write_out.write("##Original File: " + genome_gff.name + "\n##Additional Tool: "+ tool + '\n')
    for pos, data in new_ORFs.items():
        pos_ = pos.split(',')
        start = pos_[0]
        stop = pos_[-1]
        strand = data[0]
        if len(data) == 3: # Addition Tool ORFs have start and stop codons
            type = tool
            entry = (genome + '\t' + type + '\tORF\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Predicted_Additional_ORF' + '\n')
        else:
            type = 'original'
            entry = (genome + '\t' + type + '\tORF\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation' + '\n')
        write_out.write(entry)

def gff_adder(genome,tool,parameters,overlap,output_file): # Only works for single contig genome
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
    orfs_To_Remove = []

    for orf in orfs.keys():
        o_Start = int(orf.split(',')[0])
        o_Stop = int(orf.split(',')[1])
        orf_Set = set(range(int(o_Start), int(o_Stop) + 1))
        for gene in gold_standard.keys():  # Very ineffecient
            g_Start = int(gene.split(',')[0])
            g_Stop = int(gene.split(',')[1])
            gene_Set = set(range(int(g_Start), int(g_Stop) + 1))
            cov = len(orf_Set.intersection(gene_Set))
            if cov >= overlap:
                orfs_To_Remove.append(str(o_Start)+','+str(o_Stop))
            if g_Start > o_Stop:
                break
    for orf_Key in orfs_To_Remove: # Remove ORFs from out of frame if ORF was correctly matched to another Gene
        if orf_Key in orfs:
            del orfs[orf_Key]
    #########################################################
    new_ORFs = {**gold_standard, **orfs}
    new_ORFs = sortORFs(new_ORFs)
    gff_writer(genome,tool,output_file,new_ORFs,genome_gff,)

if __name__ == "__main__":
    gff_adder(**vars(args))

    print("Complete")
