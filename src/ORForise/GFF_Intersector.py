from importlib import import_module
import argparse
import collections
from datetime import date
import sys
try:
    from utils import *
except ImportError:
    from .utils import *

################################


def gff_writer(genome_ID, genome_DNA,reference_annotation, reference_tool, ref_gene_set, additional_annotation, additional_tool, genes_To_Keep, output_file):
    write_out = open(output_file, 'w')
    write_out.write("##gff-version\t3\n#\tGFF-Intersector\n#\tRun Date:" + str(date.today()) + '\n')
    write_out.write("##Genome DNA File:" + genome_DNA + '\n')
    write_out.write("##Original File: " + reference_annotation + "\n##Intersecting File: " + additional_annotation + '\n')
    for pos, data in genes_To_Keep.items():
        pos_ = pos.split(',')
        start = pos_[0]
        stop = pos_[-1]
        strand = data[0]
        type = 'original'
        entry = (
                    genome_ID + '\t' + type + '\t' + data[2] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation;Coverage=' + str(
                data[1]) + '\n')
        write_out.write(entry)


def comparator(options):  # Only works for single contig genome
    genome_seq = ""
    with open(options.genome_DNA, 'r') as genome_fasta:
        for line in genome_fasta:
            line = line.replace("\n", "")
            if not line.startswith('>'):
                genome_seq += str(line)
            else:
                genome_ID = line.split()[0].replace('>', '')
    ###########################################
    if not options.reference_tool:  # IF using Ensembl for comparison
        ref_genes = collections.OrderedDict()  # Order is important
        count = 0
        with open(options.reference_annotation, 'r') as genome_gff:
            for line in genome_gff:
                line = line.split('\t')
                try:
                    if 'CDS' in options.gene_ident and len(options.gene_ident) == 1:
                        if "CDS" in line[2] and len(line) == 9:
                            start = int(line[3])
                            stop = int(line[4])
                            strand = line[6]
                            pos = str(start) + ',' + str(stop)
                            ref_genes.update({pos: [strand, 'ref', 'CDS']})
                            count += 1
                    else:
                        gene_types = options.gene_ident.split(',')
                        if any(gene_type in line[2] for gene_type in gene_types):  # line[2] for normalrun
                            start = int(line[3])
                            stop = int(line[4])
                            strand = line[6]
                            pos = str(start) + ',' + str(stop)
                            ref_genes.update(
                                {pos: [strand, 'ref', line[2]]})  # Report what type of gene/rRNA etc we have here
                            count += 1
                except IndexError:
                    continue
    else:  # IF using a tool as reference
        try:
            reference_tool_ = import_module('Tools.' + options.reference_tool + '.' + options.reference_tool,
                                             package='my_current_pkg')
        except ModuleNotFoundError:
            try:
                reference_tool_ = import_module('ORForise.Tools.' + options.reference_tool + '.' + options.reference_tool,
                                             package='my_current_pkg')
            except ModuleNotFoundError:
                sys.exit("Tool not available")
        reference_tool_ = getattr(reference_tool_, options.reference_tool)
        ############ Reformatting tool output for ref_genes
        ref_genes = reference_tool_(options.reference_annotation, genome_seq)
    ref_gene_set = list(ref_genes.keys())
    ############################## Get Add'l
    try:
        additional_tool_ = import_module('Tools.' + options.additional_tool + '.' + options.additional_tool,
                                        package='my_current_pkg')
    except ModuleNotFoundError:
        try:
            additional_tool_ = import_module('ORForise.Tools.' + options.additional_tool + '.' + options.additional_tool,
                                            package='my_current_pkg')
        except ModuleNotFoundError:
            sys.exit("Tool not available")
    additional_tool_ = getattr(additional_tool_, options.additional_tool)
    additional_orfs = additional_tool_(options.additional_annotation, genome_seq)
    ##############################


    genes_To_Keep = collections.OrderedDict()

    if options.coverage == 100:
        for orf, data in additional_orfs.items():
            o_Start = int(orf.split(',')[0])
            o_Stop = int(orf.split(',')[1])
            o_Strand = data[0]
            try:
                if ref_genes[str(o_Start) + ',' + str(o_Stop)][2] == "CDS" : # Make sure 100% match and is also CDS
                    genes_To_Keep.update({str(o_Start) + ',' + str(o_Stop): [o_Strand, options.coverage,"CDS"]})  # o_ and g_ would be the same here
            except KeyError:
                continue
    else:
        for orf, data in additional_orfs.items():  # Currently allows ORF to be bigger than Gene
            o_Start = int(orf.split(',')[0])
            o_Stop = int(orf.split(',')[1])
            o_Strand = data[0]
            orf_Set = set(range(int(o_Start), int(o_Stop) + 1))
            for gene, g_data in ref_genes.items():  # Very ineffecient
                g_Start = int(gene.split(',')[0])
                g_Stop = int(gene.split(',')[1])
                g_Strand = g_data[0]
                gene_Set = set(range(int(g_Start), int(g_Stop) + 1))
                overlap = len(orf_Set.intersection(gene_Set))
                cov = 100 * float(overlap) / float(len(gene_Set))
                if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and cov >= options.coverage:
                    genes_To_Keep.update({str(g_Start) + ',' + str(g_Stop): [g_Strand, int(cov),g_data[2]]})
                if g_Start > o_Stop:
                    break
    #########################################################
    #### Currently, only CDSs are filtered
    for gene, g_data in ref_genes.items():  # Very ineffecient
        if "CDS" not in g_data[2]:
            g_Start = int(gene.split(',')[0])
            g_Stop = int(gene.split(',')[1])
            g_Strand = g_data[0]
            genes_To_Keep.update({str(g_Start) + ',' + str(g_Stop): [g_Strand, "N/A",g_data[2]]})
    genes_To_Keep = sortORFs(genes_To_Keep)
    gff_writer(genome_ID, options.genome_DNA, options.reference_annotation, options.reference_tool, ref_gene_set, options.additional_annotation, options.additional_tool, genes_To_Keep, options.output_file)

def main():
    print("Thank you for using ORForise\nPlease report any issues to: https://github.com/NickJD/ORForise/issues\n#####")

    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': GFF-Intersector Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-dna', dest='genome_DNA', required=True, help='Genome DNA file (.fa) which both annotations '
                                                                    'are based on')
    required.add_argument('-ref', dest='reference_annotation', required=True,
                        help='Which reference annotation file to use as reference?')
    required.add_argument('-at', dest='additional_tool', required=True,
                        help='Which format to use for additional annotation?')
    required.add_argument('-add', dest='additional_annotation', required=True,
                        help='Which annotation file to add to reference annotation?')
    required.add_argument('-o', dest='output_file', required=True,
                        help='Output filename')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-rt', dest='reference_tool', required=False,
                        help='Which tool format to use as reference? - If not provided, will default to '
                             'standard Ensembl GFF format, can be Prodigal or any of the other tools available')
    optional.add_argument('-gi', dest='gene_ident', default='CDS', required=False,
                        help='Identifier used for extraction of "genic" regions from reference annotation '
                             '"CDS,rRNA,tRNA": Default for is "CDS"')
    optional.add_argument('-cov', dest='coverage', default=100, type=int, required=False,
                        help='Percentage coverage of reference annotation needed to confirm intersection'
                             ' - Default: 100 == exact match')

    options = parser.parse_args()
    comparator(options)



if __name__ == "__main__":
    main()
    print("Complete")
