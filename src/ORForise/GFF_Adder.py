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


def gff_writer(options,genome_ID, genome_DNA, reference_annotation, reference_tool, ref_gene_set, additional_annotation, additional_tool, combined_ORFs, output_file):
    write_out = open(output_file, 'w')

    #write_out.write('##sequence-region ' + genome_ID + ' 1 ' + str(len(genome_DNA)) + '\n')
    write_out.write("##gff-version\t3\n#\tGFF-Adder\n#\tRun Date:" + str(date.today()) + '\n')
    write_out.write("##Genome DNA File:" + genome_DNA + '\n')
    write_out.write("##Original File: " + reference_annotation + "\n##Additional File: " + additional_annotation + '\n')


    #meta counts
    Ref_Only = 0
    Ref_Combined = collections.defaultdict(int)
    Non_Ref_Combined = collections.defaultdict(int)


    for pos, data in combined_ORFs.items():
        pos_ = pos.split(',')
        start = pos_[0]
        stop = pos_[-1]
        strand = data[0]
        length = int(stop) - int(start)
        additional_annotation_info = ''
        tools = additional_tool.split(',')
        matched_tools = ''
        matching = []
        matched = False
        for tool in tools:

            try:
                if options.mark_consensus == True:
                    match = [s for s in data if tool in s]
                    matching.append(match[0].replace('\n', '').replace('ID=',''))
                else:
                    match = [s for s in data if tool in s]
                    matching.append(match[0].replace('\n', '').replace('ID=',''))
                if matching:
                    matched = True

                matched_tools += tool + ','
            except Exception as e:
                if options.verbose == True:
                    print("Exception - (No matching annotation) : " + str(e))
                continue
        #temporary verbose fix
        additional_annotation_info = 'ID='
        if len(match) >1:
            for match in matching:
                additional_annotation_info += match+'|'
            additional_annotation_info = additional_annotation_info[:-1]
        elif len(match) == 1:
            additional_annotation_info += matching[0].replace('Prokka|','').replace('GeneMark_S_2|','')

        matching = None

        if pos not in ref_gene_set:  # Check if ref or additional
            type = matched_tools[:-1]
            Non_Ref_Combined[len(matched_tools.split(','))] += 1
            if options.clean == False:
                entry = (genome_ID + '\t' + type + '\t' + data[3] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Additional_Annotations;' + additional_annotation_info  +  '\n')
            else:
                entry = (genome_ID + '\t' + type + '\t' + data[3] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\t' +  additional_annotation_info + '\n')

        else:
            data[3] = data[3].replace('\n', '')#.replace('ID=', '')
            if not additional_annotation_info:
                Ref_Only += 1
                type = reference_annotation.split('/')[-1].split('.')[0]
                if options.clean == False:
                    entry = (genome_ID + '\t' + type + '\t' + data[2] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation;' + data[3] + '\n')
                else:
                    entry = (genome_ID + '\t' + type + '\t' + data[2] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\t' + data[3] + '\n')
            else:
                Ref_Combined[len(matched_tools.split(','))] +=1
                type = reference_annotation.split('/')[-1].split('.')[0]
                if options.clean == False:
                    entry = (genome_ID + '\t' + type + '\t' + data[2] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\tID=Original_Annotation;' + data[3] +
                    ';Matched_Annotations=' + additional_annotation_info + '\n')
                else:
                    entry = (genome_ID + '\t' + type + '\t' + data[2] + '\t' + start + '\t' + stop + '\t.\t' + strand + '\t.\t' + data[3] + '\n')
        write_out.write(entry)

    if options.output_meta == True:
        meta_out = open(output_file.replace('.gff','_Meta.txt'),'w')
        meta_out.write('Reference Only Genes: ' +  str(Ref_Only) + '\n')
        Ref_Combined_counter = collections.Counter(Ref_Combined)
        meta_out.write('Reference Combined Genes: ' + str(Ref_Combined_counter) + '\n')
        Non_Ref_Combined_counter = collections.Counter(Non_Ref_Combined)
        meta_out.write('Non_Reference Combined Genes: ' + str(Non_Ref_Combined_counter) + '\n')

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
                            info = line[8]
                            ref_genes.update({pos:[strand,'ref','CDS',info]})
                            count += 1
                    else:
                        gene_types = options.gene_ident.split(',')
                        if any(gene_type in line[2] for gene_type in gene_types):  # line[2] for normalrun
                            start = int(line[3])
                            stop = int(line[4])
                            strand = line[6]
                            pos = str(start) + ',' + str(stop)
                            info = line[8]
                            ref_genes.update({pos: [strand, 'ref',line[2],info]}) #Report what type of gene/rRNA etc we have here
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
    ref_genes = sortORFs(ref_genes)
    ref_gene_set = list(ref_genes.keys())
    ################ Get Additional Tool'
    # if 'StORF_Reporter' == options.additional_tool:
    #     additional_tool = 'StORF_Reporter'
    additional_annotations = collections.OrderedDict()
    tool_count = 0
    for tool in options.additional_tool.split(','):
        try:
            additional_tool_ = import_module('Tools.' + tool + '.' + tool,
                                            package='my_current_pkg')
        except ModuleNotFoundError:
            try:
                additional_tool_ = import_module('ORForise.Tools.' + tool + '.' + tool,
                                                package='my_current_pkg')
            except ModuleNotFoundError:
                sys.exit("Tool not available")
        additional_tool_ = getattr(additional_tool_, tool)
        current_additional_orfs = additional_tool_(options.additional_annotation.split(',')[tool_count], genome_seq,options.gene_ident)#,gene_ident=options.gene_ident)
        tool_count += 1
        orfs_to_remove = []
        for orf in current_additional_orfs.keys():
            o_start = int(orf.split(',')[0])
            o_stop = int(orf.split(',')[1])
            orf_set = set(range(int(o_start), int(o_stop) + 1))
            for pos, details in ref_genes.items():  # Loop through each gene to compare against predicted ORFs - Slow
                g_start = int(pos.split(',')[0])
                g_stop = int(pos.split(',')[1])

                gene_set = set(range(int(g_start), int(g_stop) + 1))
                cov = len(orf_set.intersection(gene_set))
                if g_start > o_stop:
                    break
                if cov >= options.overlap:
                    orfs_to_remove.append(str(o_start) + ',' + str(o_stop))
                    ref_genes[pos].append(current_additional_orfs[orf][4]) # record overlap
                    break
            try:
                for pos, details in additional_annotations.items():
                    a_start = int(pos.split(',')[0])
                    a_stop = int(pos.split(',')[1])
                    add_set = set(range(int(a_start), int(a_stop) + 1))
                    cov = len(orf_set.intersection(add_set))
                    if a_start > a_stop:
                        break
                    if cov >= options.overlap:
                        #orfs_to_remove.append(str(a_start) + ',' + str(a_stop))
                        additional_annotations[pos].append(current_additional_orfs[orf][4]) # record overlap
                        break
            except:
                break

        for orf_key in orfs_to_remove:  # Remove ORFs from out of frame if ORF was correctly matched to another Gene
            if orf_key in current_additional_orfs:
                del current_additional_orfs[orf_key]
        additional_annotations.update(current_additional_orfs)
    #########################################################
    combined_ORFs = {**ref_genes, **additional_annotations}
    combined_ORFs = sortORFs(combined_ORFs)

    if not options.reference_tool:
        options.reference_tool = 'Reference_Annotation'
    gff_writer(options, genome_ID, options.genome_DNA, options.reference_annotation, options.reference_tool, ref_gene_set, options.additional_annotation, options.additional_tool, combined_ORFs, options.output_file)


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
                        help='Which format to use for additional annotation? - Can provide multiple annotations (Tool1,Tool2)')
    required.add_argument('-add', dest='additional_annotation', required=True,
                        help='Which annotation file to add to reference annotation? - Can provide multiple annotations (1.GFF,2.GFF)')
    required.add_argument('-o', dest='output_file', required=True,
                        help='Output filename')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-rt', dest='reference_tool', required=False,
                        help='Which tool format to use as reference? - If not provided, will default to the '
                             'standard GFF format and will only look for "CDS" features')
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                        help='Identifier used for identifying genomic features in reference annotation "CDS,rRNA,tRNA"')
    optional.add_argument('-mc', dest='mark_consensus', default=False, type=bool, required=False,
                        help='Default - False: Mark reference annotations which where present in the additional tool annotation')
    optional.add_argument('-c', dest='clean', default=False, type=bool, required=False,
                        help='Default - False: Do not mark 9th column with "Original/Matched/Additional tag"')
    optional.add_argument('-meta', dest='output_meta', default=False, type=bool, required=False,
                        help='Default - False: Output metadata file')
    optional.add_argument('-olap', dest='overlap', default=50, type=int, required=False,
                        help='Maximum overlap between reference and additional genic regions (CDS,rRNA etc) - Default: 50 nt')

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-v', dest='verbose', default='False', type=eval, choices=[True, False],
                      help='Default - False: Print out runtime status')

    options = parser.parse_args()

    gff_adder(options)



if __name__ == "__main__":
    main()
    print("Complete")
