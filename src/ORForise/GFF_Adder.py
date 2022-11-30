from importlib import import_module
import argparse
import collections
from datetime import date
import sys
try:
    from utils import *
except ImportError:
    from .utils import *


########################################


def gff_writer(genome_ID, genome_DNA, reference_annotation, reference_tool, ref_gene_set, additional_annotation, additional_tool, combined_ORFs, output_file):
    write_out = open(output_file, 'w')
    write_out.write('##sequence-region ' + genome_ID + ' 1 ' + str(len(genome_DNA)) + '\n')
    write_out.write("##gff-version\t3\n#\tGFF-Adder\n#\tRun Date:" + str(date.today()) + '\n')
    write_out.write("##Genome DNA File:" + genome_DNA + '\n')
    write_out.write("##Original File: " + reference_annotation + "\n##Additional File: " + additional_annotation + '\n')
    storf_num = 0
    for pos, data in combined_ORFs.items():
        pos_ = pos.split(',')
        start = pos_[0]
        stop = pos_[-1]
        strand = data[0]
        length = int(stop) - int(start)
        if pos not in ref_gene_set:  # Check if ref or additional
            type = additional_tool
            entry = (
                        genome_ID + '\t' + type + '\tCDS\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Additional_Annotation_StORF_' + str(storf_num) + ';Length=' + str(length) +  '\n')
            storf_num +=1
        else:
            type = reference_tool
            entry = (
                        genome_ID + '\t' + type + '\t' + data[2] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation_' + data[1] + ';Length=' + str(length) +  '\n')
        write_out.write(entry)


def gff_adder(options):#genome_DNA, reference_tool, reference_annotation, additional_tool, additional_annotation, gene_ident, overlap, output_file):  # Only works for single contig genome
    genome_seq = ""
    with open(options.genome_DNA, 'r') as genome_fasta:
        for line in genome_fasta:
            line = line.replace("\n", "")
            if not line.startswith('>'):
                genome_seq += str(line)
            else:
                genome_ID = line.split()[0].replace('>','')
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
                            pos = str(start)+','+str(stop)
                            ref_genes.update({pos:[strand,'ref','CDS']})
                            count += 1
                    else:
                        gene_types = options.gene_ident.split(',')
                        if any(gene_type in line[2] for gene_type in gene_types):  # line[2] for normalrun
                            start = int(line[3])
                            stop = int(line[4])
                            strand = line[6]
                            pos = str(start) + ',' + str(stop)
                            ref_genes.update({pos: [strand, 'ref',line[2]]}) #Report what type of gene/rRNA etc we have here
                            count += 1
                except IndexError:
                    continue
    elif options.reference_tool: # IF using a tool as reference
        if 'StORF_Reporter' == options.reference_tool:
            reference_tool = 'StORF_Reporter'
        try:
            reference_tool_ = import_module('Tools.' + reference_tool + '.' + reference_tool,
                                             package='my_current_pkg')
        except ModuleNotFoundError:
            try:
                reference_tool_ = import_module('ORForise.Tools.' + reference_tool + '.' + reference_tool,
                                             package='my_current_pkg')
            except ModuleNotFoundError:
                sys.exit("Tool not available")
        reference_tool_ = getattr(reference_tool_, reference_tool)
        ############ Reformatting tool output for ref_genes
        ref_genes = reference_tool_(reference_annotation=options.reference_annotation, genome_seq=genome_seq,gene_ident=options.options.gene_ident)
    ref_gene_set = list(ref_genes.keys())
    ################ Get Additional Tool'
    # if 'StORF_Reporter' == options.additional_tool:
    #     additional_tool = 'StORF_Reporter'
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
    additional_orfs = additional_tool_(options.additional_annotation, genome_seq)#,gene_ident=options.gene_ident)
    orfs_to_remove = []
    for orf in additional_orfs.keys():
        o_start = int(orf.split(',')[0])
        o_stop = int(orf.split(',')[1])
        orf_set = set(range(int(o_start), int(o_stop) + 1))
        for pos, details in ref_genes.items():  # Loop through each gene to compare against predicted ORFs
            g_start = int(pos.split(',')[0])
            g_stop = int(pos.split(',')[1])
            gene_set = set(range(int(g_start), int(g_stop) + 1))
            cov = len(orf_set.intersection(gene_set))
            if cov >= options.overlap:
                orfs_to_remove.append(str(o_start) + ',' + str(o_stop))
            if g_start > o_stop:
                break
    for orf_key in orfs_to_remove:  # Remove ORFs from out of frame if ORF was correctly matched to another Gene
        if orf_key in additional_orfs:
            del additional_orfs[orf_key]
    #########################################################
    combined_ORFs = {**ref_genes, **additional_orfs}
    combined_ORFs = sortORFs(combined_ORFs)

    if not options.reference_tool:
        options.reference_tool = 'Reference_Annotation'
    gff_writer(genome_ID, options.genome_DNA, options.reference_annotation, options.reference_tool, ref_gene_set, options.additional_annotation, options.additional_tool, combined_ORFs, options.output_file)


def main():
    print("Thank you for using ORForise\nPlease report any issues to: https://github.com/NickJD/ORForise/issues\n#####")

    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': GFF-Adder Run Parameters.')
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
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                        help='Identifier used for identifying genomic features in reference annotation "CDS,rRNA,tRNA"')
    optional.add_argument('-olap', dest='overlap', default=50, type=int, required=False,
                        help='Maximum overlap between reference and additional genic regions (CDS,rRNA etc) - Default: 50 nt')

    options = parser.parse_args()

    gff_adder(options)


if __name__ == "__main__":
    main()
    print("Complete")
