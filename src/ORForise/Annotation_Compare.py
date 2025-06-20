from importlib import import_module
import argparse
import collections
import sys,os
import gzip,csv

try:
    from Comparator import tool_comparison
except ImportError:
    from .Comparator import tool_comparison

try:
    from utils import *
except ImportError:
    from ORForise.utils import *

##########################

def comparator(options):
    # with open(options.genome_DNA, mode='r') as genome:
    #     genome_Seq = "".join(line.rstrip() for line in genome if not line.startswith('>'))
    ##############################################
    # if not options.reference_tool:  # IF using Ensembl for comparison
    #     ref_genes = collections.OrderedDict()  # Order is important
    #     count = 0
    #     with open(options.reference_annotation, 'r') as genome_gff:
    #         for line in genome_gff:
    #             line = line.split('\t')
    #             try:
    #                 if "CDS" in line[2] and len(line) == 9:
    #                     start = int(line[3])
    #                     stop = int(line[4])
    #                     strand = line[6]
    #                     gene_details = [start,stop,strand]
    #                     ref_genes.update({count:gene_details})
    #                     count += 1
    #             except IndexError:
    #                 continue
    #     ref_genes = sortGenes(ref_genes) # sorted GFF refernce
    # else:  # IF using a tool as reference
    #     try:
    #         reference_tool_ = import_module('Tools.' + options.reference_tool + '.' + options.reference_tool,
    #                                          package='my_current_pkg')
    #     except ModuleNotFoundError:
    #         try:
    #             reference_tool_ = import_module('ORForise.Tools.' + options.reference_tool + '.' + options.reference_tool,
    #                                          package='my_current_pkg')
    #         except ModuleNotFoundError:
    #             sys.exit("Tool not available")
    try:
        try:  # Detect whether fasta/gff files are .gz or text and read accordingly
            fasta_in = gzip.open(options.genome_dna, 'rt')
            dna_regions = fasta_load(fasta_in)
        except:
            fasta_in = open(options.genome_dna, 'r', encoding='unicode_escape')
            dna_regions = fasta_load(fasta_in)
        try:
            gff_in = gzip.open(options.reference_annotation, 'rt')
            dna_regions = gff_load(options, gff_in, dna_regions)
        except:
            gff_in = open(options.reference_annotation, 'r', encoding='unicode_escape')
            dna_regions = gff_load(options, gff_in, dna_regions)
    except AttributeError:
        sys.exit("Attribute Error:\nStORF'ed GFF probably already exists - Must be deleted before running (-overwrite)")
    except FileNotFoundError:
        split_path = options.gff.split(os.sep)
        sys.exit("Directory '" + split_path[-2] + "' missing fna/gff files")
        ###############################################
    #print("here")
    #reference_tool_ = getattr(reference_tool_, options.reference_tool)
    ############ Reformatting tool output for ref_genes
    #ref_genes_tmp = reference_tool_(options.reference_annotation, dna_regions)
    #ref_genes = collections.OrderedDict()
    #for i, (pos, details) in enumerate(ref_genes_tmp.items()):
    #    pos = pos.split(',')
    #    ref_genes.update({i:[pos[0],pos[1],details[0]]})
    ####

    total_ref_genes = sum(
        len(v[2]) if isinstance(v[2], (list, tuple, set, dict, str)) else 1 for v in dna_regions.values())

    #############################################
    try:
        tool_ = import_module('Tools.' + options.tool + '.' + options.tool, package='my_current_pkg')
    except ModuleNotFoundError:
        try:
            tool_ = import_module('ORForise.Tools.' + options.tool + '.' + options.tool, package='my_current_pkg')
        except ModuleNotFoundError:
            sys.exit("Tool not available - Did you get the name right?")
    tool_ = getattr(tool_, options.tool)
    all_orfs = tool_(options.tool_prediction, dna_regions)
    results = tool_comparison(all_orfs, dna_regions, options.verbose)
    ############## Printing to std-out and optional csv file
    print('Genome Used: ' + str(options.genome_dna.split('/')[-1]))
    if options.reference_tool:
        print('Reference Tool Used: ' + str(options.reference_tool))
    else:
        print('Reference Used: ' + str(options.reference_annotation.split('/')[-1]))
    print('Tool Compared: ' + str(options.tool))
    print('Total Number of Reference Genes: ' + str(total_ref_genes))
    print('Number of Contigs: ' + str(len(dna_regions)))
    if options.outname:
        with open(options.outname, 'w', encoding='utf-8') as out_file:
            out_file.write('Genome Used: ' + str(options.genome_dna.split('/')[-1]) + '\n')
            if options.reference_tool:
                out_file.write('Reference Tool Used: ' + str(options.reference_tool) + '\n')
            else:
                out_file.write('Reference Used: ' + str(options.reference_annotation.split('/')[-1]) + '\n')
            out_file.write('Tool Compared: ' + str(options.tool) + '\n')
            out_file.write('Total Number of Reference Genes: ' + str(total_ref_genes) + '\n')
            out_file.write('Number of Contigs: ' + str(len(dna_regions)) + '\n')
    else:
        with open(options.outname, 'w', encoding='utf-8'):
            pass
    for dna_region, result in results.items():
        num_current_genes = len(dna_regions[dna_region][2])
        print("These are the results for: " + dna_region + '\n')
        ############################################# To get default output filename from input file details
        genome_name = options.reference_annotation.split('/')[-1].split('.')[0]
        metric_description = list(result['pred_metrics'].keys())
        metrics = list(result['pred_metrics'].values())
        rep_metric_description = list(result['rep_metrics'].keys())
        rep_metrics = list(result['rep_metrics'].values())


        print('Current Contig: ' + str(dna_region))
        print('Number of Genes: ' + str(num_current_genes))
        print('Number of ORFs: ' + str(result['pred_metrics']['Number_of_ORFs']))
        print('Perfect Matches: ' + str(result['pred_metrics']['Number_of_Perfect_Matches']) + ' [' + str(num_current_genes)+ '] - '+ format(100 * result['pred_metrics']['Number_of_Perfect_Matches']/num_current_genes,'.2f')+'%')
        print('Partial Matches: ' + str(len(result['pred_metrics']['partial_Hits'])) + ' [' + str(num_current_genes)+ '] - '+ format(100 * len(result['pred_metrics']['partial_Hits'])/num_current_genes,'.2f')+'%')
        print('Missed Genes: ' + str(len(result['rep_metrics']['genes_Undetected'])) + ' [' + str(num_current_genes)+ '] - '+ format(100 * len(result['rep_metrics']['genes_Undetected'])/num_current_genes,'.2f')+'%')
        print('Unmatched ORFs: ' + str(len(result['pred_metrics']['unmatched_ORFs'])) + ' [' + str(num_current_genes)+ '] - '+ format(100 * len(result['pred_metrics']['unmatched_ORFs'])/num_current_genes,'.2f')+'%')
        print('Multi-matched ORFs: ' + str(len(result['pred_metrics']['multi_Matched_ORFs'])) + ' [' + str(num_current_genes)+ '] - '+ format(100 * len(result['pred_metrics']['multi_Matched_ORFs'])/num_current_genes,'.2f')+'%')

        with open(options.outname, 'a', newline='\n', encoding='utf-8') as out_file:
            out_file.write('Current Contig: ' + str(dna_region) + '\n')
            out_file.write('Number of Genes: ' + str(num_current_genes) + '\n')



        if options.outname:
            with open(options.outname, 'a', newline='\n', encoding='utf-8') as out_file:
                tool_out = csv.writer(out_file, quoting=csv.QUOTE_NONE, escapechar=" ")

                # # Representative Metrics
                tool_out.writerow(['Representative_Metrics:'])
                tool_out.writerow(rep_metric_description)
                tool_out.writerow(rep_metrics)
                # All Metrics
                tool_out.writerow(['Prediction_Metrics:'])
                tool_out.writerow(list(result['pred_metrics'].keys()))
                tool_out.writerow(list(result['pred_metrics'].values()))
                # Coverage
                tool_out.writerow(['Reference_CDS_Gene_Coverage_of_Genome'])
                tool_out.writerow([result.get('gene_Coverage_Genome', 'N/A')])
                tool_out.writerow(['Predicted_CDS_Coverage_of_Genome'])
                tool_out.writerow([result.get('orf_Coverage_Genome', 'N/A')])
                tool_out.writerow(['Matched_Predicted_CDS_Coverage_of_Genome'])
                tool_out.writerow([result.get('matched_ORF_Coverage_Genome', 'N/A')])
                # Start/Stop Differences
                tool_out.writerow(['Start_Position_Difference:'])
                tool_out.writerow(result.get('start_Difference', []))
                tool_out.writerow(['Stop_Position_Difference:'])
                tool_out.writerow(result.get('stop_Difference', []))
                # Alternative Starts/Stops
                tool_out.writerow(['Alternative_Starts_Predicted:'])
                tool_out.writerow(result.get('other_Starts', []))
                tool_out.writerow(['Alternative_Stops_Predicted:'])
                tool_out.writerow(result.get('other_Stops', []))
                # Undetected Gene Metrics
                tool_out.writerow(['Undetected_Gene_Metrics:'])
                tool_out.writerow([
                    'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'
                ])
                tool_out.writerow(result.get('undetected_Gene_Metrics', []))
                # Perfect Matches
                tool_out.writerow(['Perfect_Match_Genes:'])
                for key, value in result.get('perfect_Matches', {}).items():
                    key_parts = key.split(',')
                    id = f">{genome_name}_{key_parts[0]}_{key_parts[1]}_{key_parts[2]}"
                    tool_out.writerow([id + '\n' + value + '\n'])
                # Partial Matches
                tool_out.writerow(['Partial_Match_Genes:'])
                for key, seqs in result.get('partial_Hits', {}).items():
                    key_parts = key.split(';')
                    gene_Seq = seqs[0]
                    orf_Seq = seqs[1]
                    partial = f"{key_parts[0]}\n{gene_Seq}\n{key_parts[1]}\n{orf_Seq}\n"
                    tool_out.writerow([partial])
                # Missed Genes
                tool_out.writerow(['\nMissed_Genes:'])
                for key, value in result.get('genes_Undetected', {}).items():
                    key_parts = key.split(',')
                    id = f">{genome_name}_{key_parts[0]}_{key_parts[1]}_{key_parts[2]}"
                    tool_out.writerow([id + '\n' + value + '\n'])
                # Unmatched ORF Metrics
                tool_out.writerow(['\nPredicted_CDSs_Without_Corresponding_Gene_In_Reference_Metrics:'])
                tool_out.writerow([
                    'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'
                ])
                tool_out.writerow(result.get('unmatched_ORF_Metrics', []))
                # Unmatched ORFs
                tool_out.writerow(['Predicted_CDS_Without_Corresponding_Gene_in_Reference:'])
                for key, value in result.get('unmatched_ORFs', {}).items():
                    key_parts = key.split(',')
                    id = f">{options.tool}_{key_parts[0]}_{key_parts[1]}_{key_parts[2]}"
                    tool_out.writerow([id + '\n' + value])
                # Multi-matched ORFs
                tool_out.writerow(['\nPredicted_CDSs_Which_Detected_more_than_one_Gene:'])
                try:
                    for key, value in result.get('multi_Matched_ORFs', {}).items():
                        key_parts = key.split(',')
                        multi = f"Predicted_CDS:{key_parts[0]}-{key_parts[1]}_Genes:{'|'.join(value)}"
                        tool_out.writerow([multi])
                except Exception:
                    pass
        #
        # if options.outname:
        #     with open(options.outname, 'wa', newline='\n', encoding='utf-8') as out_file:  # Clear write out of report
        #         tool_out = csv.writer(out_file, quoting=csv.QUOTE_NONE, escapechar=" ")
        #         tool_out.writerow(['Representative_Metrics:'])
        #         tool_out.writerow(rep_metric_description)
        #         tool_out.writerow(rep_metrics)
        #         tool_out.writerow(['All_Metrics:'])
        #         tool_out.writerow(metric_description)
        #         tool_out.writerow(metrics)
        #         tool_out.writerow(['Reference_CDS_Gene_Coverage_of_Genome'])
        #         tool_out.writerow([result['gene_Coverage_Genome']])
        #         tool_out.writerow(['Predicted_CDS_Coverage_of_Genome'])
        #         tool_out.writerow([result['orf_Coverage_Genome']])
        #         tool_out.writerow(['Matched_Predicted_CDS_Coverage_of_Genome'])
        #         tool_out.writerow([result['matched_ORF_Coverage_Genome']])
        #         tool_out.writerow(['Start_Position_Difference:'])
        #         tool_out.writerow(result['start_Difference'])
        #         tool_out.writerow(['Stop_Position_Difference:'])
        #         tool_out.writerow(result['stop_Difference'])
        #         tool_out.writerow(['Alternative_Starts_Predicted:'])
        #         tool_out.writerow(result['other_Starts'])
        #         tool_out.writerow(['Alternative_Stops_Predicted:'])
        #         tool_out.writerow(result['other_Stops'])
        #         tool_out.writerow(['Undetected_Gene_Metrics:'])
        #         tool_out.writerow([
        #                               'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'])
        #         tool_out.writerow(result['undetected_Gene_Metrics'])
        #         ####
        #         tool_out.writerow(['Perfect_Match_Genes:'])
        #         for key, value in result['perfect_Matches'].items():
        #             key = key.split(',')
        #             id = ('>' + genome_name + '_' + key[0] + '_' + key[1] + '_' + key[2])
        #             tool_out.writerow([id + '\n' + value + '\n'])
        #         ####
        #         tool_out.writerow(['Partial_Match_Genes:'])
        #         for key, seqs in result['partial_Hits'].items():
        #             key = key.split(';')
        #             gene_Seq = seqs[0]
        #             orf_Seq = seqs[1]
        #             partial = (key[0] + '\n' + gene_Seq + '\n' + key[1] + '\n' + orf_Seq + '\n')
        #             tool_out.writerow([partial])
        #         ####
        #         tool_out.writerow(['\nMissed_Genes:'])
        #         for key, value in result['genes_Undetected'].items():
        #             key = key.split(',')
        #             id = ('>' + genome_name + '_' + key[0] + '_' + key[1] + '_' + key[2])
        #             tool_out.writerow([id + '\n' + value + '\n'])
        #         tool_out.writerow(['\nPredicted_CDSs_Without_Corresponding_Gene_In_Reference_Metrics:'])
        #         tool_out.writerow([
        #                               'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'])
        #         tool_out.writerow(result['unmatched_ORF_Metrics'])
        #         tool_out.writerow(['Predicted_CDS_Without_Corresponding_Gene_in_Reference:'])
        #         for key, value in result['unmatched_ORFs'].items():
        #             key = key.split(',')
        #             id = ('>' + options.tool + '_' + key[0] + '_' + key[1] + '_' + key[2])
        #             tool_out.writerow([id + '\n' + value])
        #         tool_out.writerow(['\nPredicted_CDSs_Which_Detected_more_than_one_Gene:'])
        #
        #         try:
        #             for key, value in result['multi_Matched_ORFs'].items():
        #                 key = key.split(',')
        #                 multi = ('Predicted_CDS:' + key[0] + '-' + key[1] + '_Genes:' + '|'.join(value))
        #                 tool_out.writerow([multi])
        #         except IndexError:
        #             pass



def main():
    print("Thank you for using ORForise\nPlease report any issues to: https://github.com/NickJD/ORForise/issues\n#####")

    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': Annotatione-Compare Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-dna', dest='genome_dna', required=True, help='Genome DNA file (.fa) which both annotations '
                                                                    'are based on')
    required.add_argument('-ref', dest='reference_annotation', required=True,
                        help='Which reference annotation file to use as reference?')
    required.add_argument('-t', dest='tool', required=True, help='Which tool to analyse? (Prodigal)')
    required.add_argument('-tp', dest='tool_prediction', required=True,
                        help='Tool genome prediction file (.gff) - Different Tool Parameters'
                             ' are compared individually via separate files')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                          help='What features to consider as genes? - Default: CDS - '
                               'Provide comma separated list of features to consider as genes (e.g. CDS,exon)')
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
