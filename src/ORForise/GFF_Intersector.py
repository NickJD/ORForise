from importlib import import_module

import argparse
import collections
from datetime import date

parser = argparse.ArgumentParser()
parser.add_argument('-g', action='store', dest='genome', required=True,
                    help='Which genome to use as Gold Standard?')
parser.add_argument('-t', action='store', dest='tool', required=True,
                    help='While tool to add to Gold Standard Annotation?')
parser.add_argument('-p', action='store', dest='parameters', required=False,
                    help='Optional parameters for prediction tool.')
parser.add_argument('-cov', action='store', dest='coverage', default=100, type=int, required=False,
                    help='ORF coverage of Gene in percentage - Default of 100 means exact match')
parser.add_argument('-o', action='store', dest='output_file', required=True,
                    help='output filename')
args = parser.parse_args()


def gff_writer(genome, tool, outfile, genes_To_Keep, genome_gff):
    write_out = open(outfile, 'w')
    write_out.write("##gff-version\t3\n#\tGFF Intersector\n#\tRun Date:" + str(date.today()) + '\n')
    write_out.write("##Original File: " + genome_gff.name + "\n##Intersecting Tool: " + tool + '\n')
    for pos, data in genes_To_Keep.items():
        pos_ = pos.split(',')
        start = pos_[0]
        stop = pos_[-1]
        strand = data[0]
        type = 'original'
        entry = (
                    genome + '\t' + type + '\tORF\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation;Coverage=' + str(
                data[1]) + '\n')
        write_out.write(entry)


def comparator(genome, tool, parameters, coverage, output_file):  # Only works for single contig genome
    genome_Seq = ""
    with open('Genomes/' + genome + '.fa', 'r') as genome_fasta:
        for line in genome_fasta:
            line = line.replace("\n", "")
            if not line.startswith('>'):
                genome_Seq += str(line)
    genome_Size = len(genome_Seq)
    ###########################################
    gold_standard = collections.OrderedDict()  # Order is important
    with open('Genomes/' + genome + '.gff', 'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "CDS" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]
                    pos = str(start) + ',' + str(stop)
                    gold_standard.update({pos: strand})
            except IndexError:
                pass

    tool_predictions = import_module('Tools.' + tool + '.' + tool)
    tool_predictions = getattr(tool_predictions, tool)
    orfs = tool_predictions(genome, parameters, genome_Seq)
    genes_To_Keep = collections.OrderedDict()

    if coverage == 100:
        for orf, data in orfs.items():
            o_Start = int(orf.split(',')[0])
            o_Stop = int(orf.split(',')[1])
            o_Strand = data[0]
            try:
                if gold_standard[str(o_Start) + ',' + str(o_Stop)]:
                    genes_To_Keep.update(
                        {str(o_Start) + ',' + str(o_Stop): [o_Strand, coverage]})  # o_ and g_ would be the same here
            except KeyError:
                continue
    else:
        for orf, data in orfs.items():  # Currently allows ORF to be bigger than Gene
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
                    genes_To_Keep.update({str(g_Start) + ',' + str(g_Stop): [g_Strand, int(cov)]})
                if g_Start > o_Stop:
                    break
    #########################################################
    gff_writer(genome, tool, output_file, genes_To_Keep, genome_gff)


if __name__ == "__main__":
    comparator(**vars(args))

    print("Complete")
