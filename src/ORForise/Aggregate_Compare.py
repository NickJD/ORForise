from importlib import import_module
import argparse
import collections
import csv
import sys

try:
    from Comparator import tool_comparison
    from utils import *
except ImportError:
    from .Comparator import tool_comparison
    from .utils import *

############################################

def comparator(options):
    genome_seq = ""
    with open(options.genome_DNA, 'r') as genome:
        for line in genome:
            line = line.replace("\n", "")
            if not line.startswith('>'):
                genome_seq += str(line)
            else:
                genome_ID = line.split()[0].replace('>','')
    ##############################################
    if not options.reference_tool:  # IF using Ensembl for comparison
        ref_genes = collections.OrderedDict()  # Order is important
        count = 0
        with open(options.reference_annotation, 'r') as genome_gff:
            for line in genome_gff:
                line = line.split('\t')
                try:
                    if "CDS" in line[2] and len(line) == 9:
                        start = int(line[3])
                        stop = int(line[4])
                        strand = line[6]
                        gene_details = [start, stop, strand]
                        ref_genes.update({count: gene_details})
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
        ref_genes_tmp = reference_tool_(options.reference_annotation, genome_seq)
        ref_genes = collections.OrderedDict()
        for i, (pos, details) in enumerate(ref_genes_tmp.items()):
            pos = pos.split(',')
            ref_genes.update({i: [pos[0], pos[1], details[0]]})
    #############################################
    # Currently only one model type can be used. (--parameters)
    aggregate_Predictions = collections.OrderedDict()
    aggregate_Tools = options.tools.split(',')
    for i, (tool) in enumerate(aggregate_Tools):
        tool_prediction = options.tool_predictions.split(',')[i]
        print(tool)
        try:
            tool_ = import_module('Tools.' + tool + '.' + tool, package='my_current_pkg')
        except ModuleNotFoundError:
            try:
                tool_ = import_module('ORForise.Tools.' + tool + '.' + tool, package='my_current_pkg')
            except ModuleNotFoundError:
                sys.exit("Tool not available")
        tool_ = getattr(tool_, tool)
        orfs = tool_(tool_prediction, genome_seq)
        aggregate_Predictions.update(orfs)

    aggregate_Predictions = sortORFs(aggregate_Predictions)
    all_Metrics, all_rep_Metrics, start_precision, stop_precision, other_starts, other_stops, perfect_Matches, missed_genes, unmatched_orfs, undetected_gene_metrics, unmatched_orf_metrics, orf_Coverage_Genome, matched_ORF_Coverage_Genome, gene_coverage_genome, multi_Matched_ORFs, partial_Hits = tool_comparison(
        ref_genes, aggregate_Predictions, genome_seq, options.verbose)
    ############################################# To get default output filename from input file details
    genome_name = options.reference_annotation.split('/')[-1].split('.')[0]
    metric_description = list(all_Metrics.keys())
    metrics = list(all_Metrics.values())
    rep_metric_description = list(all_rep_Metrics.keys())
    rep_metrics = list(all_rep_Metrics.values())
    #################################
    print('Genome Used: ' + str(options.reference_annotation.split('/')[-1]))
    if options.reference_tool:
        print('Reference Tool Used: ' + str(options.reference_tool))
    else:
        print('Reference Used: ' + str(options.reference_annotation))
    print('Tools Compared: ' + str(options.tools))
    print('Perfect Matches:' + str(len(perfect_Matches)) + '[' + str(len(ref_genes))+ '] -'+ format(100 * len(perfect_Matches)/len(ref_genes),'.2f')+'%')
    print('Partial Matches:' + str(len(partial_Hits)) + '[' + str(len(ref_genes))+ '] - '+ format(100 * len(partial_Hits)/len(ref_genes),'.2f')+'%')
    print('Missed Genes:' + str(len(missed_genes)) + '[' + str(len(ref_genes))+ '] - '+ format(100 * len(missed_genes)/len(ref_genes),'.2f')+'%')
    if options.outname:
        with open(options.outname, 'w', newline='\n',
                  encoding='utf-8') as out_file:  # Clear write out of report
            tool_out = csv.writer(out_file, quoting=csv.QUOTE_NONE, escapechar=" ")
            tool_out.writerow(['Representative_Metrics:'])
            tool_out.writerow(rep_metric_description)
            tool_out.writerow(rep_metrics)
            tool_out.writerow(['All_Metrics:'])
            tool_out.writerow(metric_description)
            tool_out.writerow(metrics)
            tool_out.writerow(['Reference_CDS_Gene_Coverage_of_Genome'])
            tool_out.writerow([gene_coverage_genome])
            tool_out.writerow(['Predicted_CDS_Coverage_of_Genome'])
            tool_out.writerow([orf_Coverage_Genome])
            tool_out.writerow(['Matched_Predicted_CDS_Coverage_of_Genome'])
            tool_out.writerow([matched_ORF_Coverage_Genome])
            tool_out.writerow(['Start_Position_Difference:'])
            tool_out.writerow(start_precision)
            tool_out.writerow(['Stop_Position_Difference:'])
            tool_out.writerow(stop_precision)
            tool_out.writerow(['Alternative_Starts_Predicted:'])
            tool_out.writerow(other_starts)
            tool_out.writerow(['Alternative_Stops_Predicted:'])
            tool_out.writerow(other_stops)
            tool_out.writerow(['Undetected_Gene_Metrics:'])
            tool_out.writerow([
                                  'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'])
            tool_out.writerow(undetected_gene_metrics)
            tool_out.writerow(['Perfect_Match_Genes:'])
            for key, value in perfect_Matches.items():
                key = key.split(',')
                id = ('>' + genome_name + '_' + key[0] + '_' + key[1] + '_' + key[2])
                tool_out.writerow([id + '\n' + value + '\n'])
            ####
            tool_out.writerow(['Partial_Match_Genes:'])
            for key, seqs in partial_Hits.items():
                key = key.split(';')
                gene_Seq = seqs[0]
                orf_Seq = seqs[1]
                partial = (key[0] + '\n' + gene_Seq + '\n' + key[1] + '\n' + orf_Seq + '\n')
                tool_out.writerow([partial])
            ####
            tool_out.writerow(['\nMissed_Genes:'])
            for key, value in missed_genes.items():
                key = key.split(',')
                id = ('>' + genome_name + '_' + key[0] + '_' + key[1] + '_' + key[2])
                tool_out.writerow([id + '\n' + value + '\n'])
            tool_out.writerow(['\nPredicted_CDSs_Without_Corresponding_Gene_In_Reference_Metrics:'])
            tool_out.writerow([
                'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'])
            tool_out.writerow(unmatched_orf_metrics)
            tool_out.writerow(['Predicted_CDS_Without_Corresponding_Gene_in_Reference:'])
            for key, value in unmatched_orfs.items():
                key = key.split(',')
                id = ('>' + tool + '_' + key[0] + '_' + key[1] + '_' + key[2])
                tool_out.writerow([id + '\n' + value])
            tool_out.writerow(['\nPredicted_CDSs_Which_Detected_more_than_one_Gene:'])

            try:
                for key, value in multi_Matched_ORFs.items():
                    key = key.split(',')
                    multi = ('Predicted_CDS:' + key[0] + '-' + key[1] + '_Genes:' + '|'.join(value))
                    tool_out.writerow([multi])
            except IndexError:
                pass


def main():
    print("Thank you for using ORForise\nPlease report any issues to: https://github.com/NickJD/ORForise/issues\n#####")

    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': Aggregate-Compare Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')

    required.add_argument('-dna', dest='genome_DNA', required=True, help='Genome DNA file (.fa) which both annotations '
                                                                    'are based on')
    required.add_argument('-t', dest='tools', required=True, help='Which tools to analyse? (Prodigal,GeneMarkS)')
    required.add_argument('-tp', dest='tool_predictions', required=True, help='Tool genome prediction file (.gff) - Provide'
                                                                         'file locations for each tool comma separated')
    required.add_argument('-ref', dest='reference_annotation', required=True,
                        help='Which reference annotation file to use as reference?')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-rt', dest='reference_tool', required=False,
                        help='What type of Annotation to compare to? -- Leave blank for Ensembl reference'
                             '- Provide tool name to compare output from two tools')

    output = parser.add_argument_group('Output')
    output.add_argument('-o', dest='outname', required=False,
                        help='Define full output filename (format is CSV) - If not provided, summary will be printed to std-out')

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-v', dest='verbose', default='False', type=eval, choices=[True, False],
                        help='Default - False: Print out runtime status')
    options = parser.parse_args()
    comparator(options)

if __name__ == "__main__":
    main()
    print("Complete")

