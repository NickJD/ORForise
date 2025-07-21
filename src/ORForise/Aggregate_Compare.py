from importlib import import_module
import argparse
import csv, os, gzip, sys


try:
    from Comparator import tool_comparison
    from utils import *
except ImportError:
    from .Comparator import tool_comparison
    from .utils import *

############################################

def comparator(options):
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
    total_ref_genes = sum(
        len(v[2]) if isinstance(v[2], (list, tuple, set, dict, str)) else 1 for v in dna_regions.values())
    #############################################
    # Collect predictions from tools
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
        ##
        orfs = tool_(tool_prediction, dna_regions)
        for current_contig in orfs:
            if current_contig not in aggregate_Predictions:
                aggregate_Predictions[current_contig] = {}
            current_orfs = orfs[current_contig]
            for key, value in current_orfs.items():
                if key in aggregate_Predictions[current_contig]:
                    aggregate_Predictions[current_contig][key][-1] += '|' + tool
                else:
                    aggregate_Predictions[current_contig][key] = value

    aggregate_ORFs = {k: sortORFs(v) for k, v in aggregate_Predictions.items()}
    results = tool_comparison(aggregate_ORFs, dna_regions, options.verbose)
    ############## Printing to std-out and optional csv file
    # Ensure the output directory exists
    os.makedirs(options.outdir, exist_ok=True)
    # Use outname as a directory, basename for files is output-outname
    base_out = os.path.join(options.outdir, f"{os.path.basename(options.outname)}")

    # Prepare to collect summary stats for all contigs
    contig_summaries = []
    ############################################# To get default output filename from input file details
    if options.outdir:
        # Ensure the output directory exists
        os.makedirs(options.outdir, exist_ok=True)
        # Use outname as a directory, basename for files is output-outname
        base_out = os.path.join(options.outdir, f"{os.path.basename(options.outname)}")
        with open(f"{base_out}_summary.txt", 'w', encoding='utf-8') as out_file:
            out_file.write('Genome Used: ' + str(options.genome_dna.split('/')[-1]) + '\n')
            if options.reference_tool:
                out_file.write('Reference Tool Used: ' + str(options.reference_tool) + '\n')
            else:
                out_file.write('Reference Used: ' + str(options.reference_annotation.split('/')[-1]) + '\n')
            out_file.write('Tool Compared: ' + str(options.tools) + '\n')
            out_file.write('Total Number of Reference Genes: ' + str(total_ref_genes) + '\n')
            out_file.write('Number of Contigs: ' + str(len(dna_regions)) + '\n')
            out_file.write(
                'Contig\tGenes\tORFs\tPerfect_Matches\tPartial_Matches\tMissed_Genes\tUnmatched_ORFs\tMulti_Matched_ORFs\n')

    for dna_region, result in results.items():
        num_current_genes = len(dna_regions[dna_region][2])
        num_orfs = result['pred_metrics']['Number_of_ORFs']
        num_perfect = result['pred_metrics']['Number_of_Perfect_Matches']
        num_partial = len(result['pred_metrics']['partial_Hits'])
        num_missed = len(result['rep_metrics']['genes_Undetected'])
        num_unmatched = len(result['pred_metrics']['unmatched_ORFs'])
        num_multi = len(result['pred_metrics']['multi_Matched_ORFs'])

        ####
        # Tool-specific stats
        tool_stats = {}
        for tool in options.tools.split(','):
            tool_stats[tool] = {
                'perfect': 0,
                'partial': 0,
                'unmatched': 0,
                'multi': 0
            }
        # Count perfect matches per tool
        for key in result['pred_metrics'].get('perfect_Matches', {}):
            for tool in options.tools.split(','):
                if tool in key:
                    tool_stats[tool]['perfect'] += 1
        # Count partial matches per tool
        for key in result['pred_metrics'].get('partial_Hits', {}):
            for tool in options.tools.split(','):
                if tool in key:
                    tool_stats[tool]['partial'] += 1
        # Count unmatched ORFs per tool
        for key in result['pred_metrics'].get('unmatched_ORFs', {}):
            for tool in options.tools.split(','):
                if tool in key:
                    tool_stats[tool]['unmatched'] += 1
        # Count multi-matched ORFs per tool
        for key in result['pred_metrics'].get('multi_Matched_ORFs', {}):
            for tool in options.tools.split(','):
                if tool in key:
                    tool_stats[tool]['multi'] += 1
        ####

        # Collect summary for this contig
        if options.outdir:
            contig_summaries.append([
                dna_region, num_current_genes, num_orfs, num_perfect, num_partial, num_missed, num_unmatched, num_multi
            ])
        ###
        num_current_genes = len(dna_regions[dna_region][2])
        print("These are the results for: " + dna_region + '\n')
        ############################################# To get default output filename from input file details
        genome_name = options.reference_annotation.split('/')[-1].split('.')[0]
        rep_metric_description, rep_metrics = get_rep_metrics(result)
        all_metric_description, all_metrics = get_all_metrics(result)

        print('Current Contig: ' + str(dna_region))
        print('Number of Genes: ' + str(num_current_genes))
        print('Number of ORFs: ' + str(result['pred_metrics']['Number_of_ORFs']))
        print('Perfect Matches: ' + str(result['pred_metrics']['Number_of_Perfect_Matches']) + ' [' + str(num_current_genes)+ '] - '+ format(100 * result['pred_metrics']['Number_of_Perfect_Matches']/num_current_genes,'.2f')+'%')
        print('Partial Matches: ' + str(len(result['pred_metrics']['partial_Hits'])) + ' [' + str(num_current_genes)+ '] - '+ format(100 * len(result['pred_metrics']['partial_Hits'])/num_current_genes,'.2f')+'%')
        print('Missed Genes: ' + str(len(result['rep_metrics']['genes_Undetected'])) + ' [' + str(num_current_genes)+ '] - '+ format(100 * len(result['rep_metrics']['genes_Undetected'])/num_current_genes,'.2f')+'%')
        print('Unmatched ORFs: ' + str(len(result['pred_metrics']['unmatched_ORFs'])) + ' [' + str(num_current_genes)+ '] - '+ format(100 * len(result['pred_metrics']['unmatched_ORFs'])/num_current_genes,'.2f')+'%')
        print('Multi-matched ORFs: ' + str(len(result['pred_metrics']['multi_Matched_ORFs'])) + ' [' + str(num_current_genes)+ '] - '+ format(100 * len(result['pred_metrics']['multi_Matched_ORFs'])/num_current_genes,'.2f')+'%')
        print('Tool breakdown:')
        for tool, stats in tool_stats.items():
            print(
                f"  {tool}: Perfect={stats['perfect']}, Partial={stats['partial']}, Unmatched={stats['unmatched']}, Multi-matched={stats['multi']}")

        if options.outdir:
            # Prepare output directory and file names for each contig
            contig_save = dna_region.replace('/', '_').replace('\\', '_')
            contig_dir = os.path.join(options.outdir, contig_save)
            os.makedirs(contig_dir, exist_ok=True)
            summary_file = os.path.join(contig_dir, "summary.txt")
            csv_file = os.path.join(contig_dir, "metrics.csv")
            perfect_fasta = os.path.join(contig_dir, "perfect_matches.fasta")
            partial_fasta = os.path.join(contig_dir, "partial_matches.fasta")
            missed_fasta = os.path.join(contig_dir, "missed_genes.fasta")
            unmatched_fasta = os.path.join(contig_dir, "unmatched_orfs.fasta")
            multi_fasta = os.path.join(contig_dir, "multi_matched_orfs.fasta")

            # Write summary to text file
            with open(summary_file, 'w', encoding='utf-8') as sf:
                sf.write('Current Contig: ' + str(dna_region) + '\n')
                sf.write('Number of Genes: ' + str(num_current_genes) + '\n')
                sf.write('Number of ORFs: ' + str(result['pred_metrics']['Number_of_ORFs']) + '\n')
                sf.write('Perfect Matches: ' + str(result['pred_metrics']['Number_of_Perfect_Matches']) + ' [' + str(
                    num_current_genes) + '] - ' + format(
                    100 * result['pred_metrics']['Number_of_Perfect_Matches'] / num_current_genes, '.2f') + '%\n')
                sf.write('Partial Matches: ' + str(len(result['pred_metrics']['partial_Hits'])) + ' [' + str(
                    num_current_genes) + '] - ' + format(
                    100 * len(result['pred_metrics']['partial_Hits']) / num_current_genes, '.2f') + '%\n')
                sf.write('Missed Genes: ' + str(len(result['rep_metrics']['genes_Undetected'])) + ' [' + str(
                    num_current_genes) + '] - ' + format(
                    100 * len(result['rep_metrics']['genes_Undetected']) / num_current_genes, '.2f') + '%\n')
                sf.write('Unmatched ORFs: ' + str(len(result['pred_metrics']['unmatched_ORFs'])) + ' [' + str(
                    num_current_genes) + '] - ' + format(
                    100 * len(result['pred_metrics']['unmatched_ORFs']) / num_current_genes, '.2f') + '%\n')
                sf.write('Multi-matched ORFs: ' + str(len(result['pred_metrics']['multi_Matched_ORFs'])) + ' [' + str(
                    num_current_genes) + '] - ' + format(
                    100 * len(result['pred_metrics']['multi_Matched_ORFs']) / num_current_genes, '.2f') + '%\n')
                sf.write('Tool breakdown:\n')
                for tool, stats in tool_stats.items():
                    sf.write(
                        f"  {tool}: Perfect={stats['perfect']}, Partial={stats['partial']}, Unmatched={stats['unmatched']}, Multi-matched={stats['multi']}\n")

            # Write metrics to CSV
            with open(csv_file, 'w', newline='\n', encoding='utf-8') as out_file:
                tool_out = csv.writer(out_file, quoting=csv.QUOTE_NONE, escapechar=" ")
                tool_out.writerow(['Representative_Metrics:'])
                tool_out.writerow(rep_metric_description.split(','))
                tool_out.writerow([*rep_metrics])
                tool_out.writerow(['Prediction_Metrics:'])
                tool_out.writerow(all_metric_description.split(','))
                tool_out.writerow([*all_metrics])
                tool_out.writerow(['Reference_CDS_Gene_Coverage_of_Genome'])
                tool_out.writerow([''.join(map(str, result['rep_metrics']['gene_Coverage_Genome']))])
                tool_out.writerow(['Predicted_CDS_Coverage_of_Genome'])
                tool_out.writerow([''.join(map(str, result['pred_metrics']['orf_Coverage_Genome']))])
                tool_out.writerow(['Matched_Predicted_CDS_Coverage_of_Genome'])
                tool_out.writerow([''.join(map(str, result['pred_metrics']['matched_ORF_Coverage_Genome']))])
                # tool_out.writerow(['Start_Position_Difference:'])
                # tool_out.writerow(result.get('start_Difference', []))
                # tool_out.writerow(['Stop_Position_Difference:'])
                # tool_out.writerow(result.get('stop_Difference', []))
                # tool_out.writerow(['Alternative_Starts_Predicted:'])
                # tool_out.writerow(result.get('other_Starts', []))
                # tool_out.writerow(['Alternative_Stops_Predicted:'])
                # tool_out.writerow(result.get('other_Stops', []))
                # tool_out.writerow(['Undetected_Gene_Metrics:'])
                # tool_out.writerow([
                #     'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'
                # ])
                # tool_out.writerow(result.get('undetected_Gene_Metrics', []))
                # tool_out.writerow(['\nPredicted_CDSs_Without_Corresponding_Gene_In_Reference_Metrics:'])
                # tool_out.writerow([
                #     'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'
                # ])
                # tool_out.writerow(result.get('unmatched_ORF_Metrics', []))

            # Write perfect matches to FASTA
            with open(perfect_fasta, 'w', encoding='utf-8') as f:
                for key, value in result['pred_metrics'].get('perfect_Matches', {}).items():
                    key_parts = key.split(',')
                    id = f">{genome_name}_{key_parts[0]}_{key_parts[1]}_{key_parts[2]}_{key_parts[5]}"
                    f.write(f"{id}\n{value}\n")

            # Write partial matches to FASTA
            with open(partial_fasta, 'w', encoding='utf-    8') as f:
                for key, value in result['pred_metrics'].get('partial_Hits', {}).items():
                    key_parts = key.split(';')
                    gene_Seq = value[0]
                    orf_Seq = value[1]
                    f.write(f">{key_parts[0]}_gene\n{gene_Seq}\n>{key_parts[1]}_orf\n{orf_Seq}\n")

            # Write missed genes to FASTA
            with open(missed_fasta, 'w', encoding='utf-8') as f:
                for key, value in result['rep_metrics'].get('genes_Undetected', {}).items():
                    key_parts = key.split(',')
                    id = f">{genome_name}_{key_parts[0]}_{key_parts[1]}_{key_parts[2]}"
                    f.write(f"{id}\n{value}\n")

            # Write unmatched ORFs to FASTA
            with open(unmatched_fasta, 'w', encoding='utf-8') as f:
                for key, value in result['pred_metrics'].get('unmatched_ORFs', {}).items():
                    key_parts = key.split(',')
                    id = f">{options.tools}_{key_parts[0]}_{key_parts[1]}_{key_parts[2]}"
                    f.write(f"{id}\n{value}\n")

            # Write multi-matched ORFs to FASTA
            with open(multi_fasta, 'w', encoding='utf-8') as f:
                for key, value in result['pred_metrics'].get('multi_Matched_ORFs', {}).items():
                    key_parts = key.split(',')
                    multi = f">Predicted_CDS:{key_parts[0]}-{key_parts[1]}_Genes:{'|'.join(value)}"
                    f.write(f"{multi}\n")

    # After all contigs, append the summary table to the main summary file
    if options.outdir and contig_summaries:
        with open(f"{base_out}_summary.txt", 'a', encoding='utf-8') as out_file:
            for row in contig_summaries:
                out_file.write('\t'.join(map(str, row)) + '\n')
            # Optionally, add overall totals
            total_genes = sum(row[1] for row in contig_summaries)
            total_orfs = sum(row[2] for row in contig_summaries)
            total_perfect = sum(row[3] for row in contig_summaries)
            total_partial = sum(row[4] for row in contig_summaries)
            total_missed = sum(row[5] for row in contig_summaries)
            total_unmatched = sum(row[6] for row in contig_summaries)
            total_multi = sum(row[7] for row in contig_summaries)
            out_file.write('\nOverall Summary:\n')
            out_file.write(f'Number of Genes: {total_genes}\n')
            out_file.write(f'Number of ORFs: {total_orfs}\n')
            out_file.write(
                f'Perfect Matches: {total_perfect} [{total_genes}] - {format(100 * total_perfect / total_genes, ".2f")}%\n')
            out_file.write(
                f'Partial Matches: {total_partial} [{total_genes}] - {format(100 * total_partial / total_genes, ".2f")}%\n')
            out_file.write(
                f'Missed Genes: {total_missed} [{total_genes}] - {format(100 * total_missed / total_genes, ".2f")}%\n')
            out_file.write(
                f'Unmatched ORFs: {total_unmatched} [{total_genes}] - {format(100 * total_unmatched / total_genes, ".2f")}%\n')
            out_file.write(
                f'Multi-matched ORFs: {total_multi} [{total_genes}] - {format(100 * total_multi / total_genes, ".2f")}%\n')

            # Calculate combined tool stats - could be optimised further
            combined_tool_stats = {tool: {'perfect': 0, 'partial': 0, 'unmatched': 0, 'multi': 0} for tool in
                                   options.tools.split(',')}
            for dna_region, result in results.items():
                for tool in options.tools.split(','):
                    # perfect
                    for key in result['pred_metrics'].get('perfect_Matches', {}):
                        if tool in key:
                            combined_tool_stats[tool]['perfect'] += 1
                    # partial
                    for key in result['pred_metrics'].get('partial_Hits', {}):
                        if tool in key:
                            combined_tool_stats[tool]['partial'] += 1
                    # unmatched
                    for key in result['pred_metrics'].get('unmatched_ORFs', {}):
                        if tool in key:
                            combined_tool_stats[tool]['unmatched'] += 1
                    # multi
                    for key in result['pred_metrics'].get('multi_Matched_ORFs', {}):
                        if tool in key:
                            combined_tool_stats[tool]['multi'] += 1
            for tool, stats in combined_tool_stats.items():
                out_file.write('\n'+
                    f"  {tool}: Perfect={stats['perfect']}, Partial={stats['partial']}, Unmatched={stats['unmatched']}, Multi-matched={stats['multi']}\n"
                )

            # Print combined metrics to stdout
            print("\nCombined metrics for all contigs:")
            print(f'Number of Genes: {total_genes}')
            print(f'Number of ORFs: {total_orfs}')
            print(
                f'Perfect Matches: {total_perfect} [{total_genes}] - {format(100 * total_perfect / total_genes, ".2f")}%')
            print(
                f'Partial Matches: {total_partial} [{total_genes}] - {format(100 * total_partial / total_genes, ".2f")}%')
            print(f'Missed Genes: {total_missed} [{total_genes}] - {format(100 * total_missed / total_genes, ".2f")}%')
            print(
                f'Unmatched ORFs: {total_unmatched} [{total_genes}] - {format(100 * total_unmatched / total_genes, ".2f")}%')
            print(
                f'Multi-matched ORFs: {total_multi} [{total_genes}] - {format(100 * total_multi / total_genes, ".2f")}%')

            print('Tool breakdown (combined):')
            for tool, stats in combined_tool_stats.items():
                print('\n'+
                    f"  {tool}: Perfect={stats['perfect']}, Partial={stats['partial']}, Unmatched={stats['unmatched']}, Multi-matched={stats['multi']}"
                )


def main():
    print("Thank you for using ORForise\nPlease report any issues to: https://github.com/NickJD/ORForise/issues\n"
          "Please Cite: https://doi.org/10.1093/bioinformatics/btab827\n"
          "#####")

    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': Aggregate-Compare Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')

    required.add_argument('-dna', dest='genome_dna', required=True, help='Genome DNA file (.fa) which both annotations '
                                                                    'are based on')
    required.add_argument('-t', dest='tools', required=True, help='Which tools to analyse? (Prodigal,GeneMarkS)')
    required.add_argument('-tp', dest='tool_predictions', required=True, help='Tool genome prediction file (.gff) - Provide'
                                                                         'file locations for each tool comma separated')
    required.add_argument('-ref', dest='reference_annotation', required=True,
                        help='Which reference annotation file to use as reference?')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                          help='What features to consider as genes? - Default: CDS - '
                               'Provide comma separated list of features to consider as genes (e.g. CDS,exon)')
    optional.add_argument('-rt', dest='reference_tool', required=False,
                        help='What type of Annotation to compare to? -- Leave blank for Ensembl reference'
                             '- Provide tool name to compare output from two tools')

    output = parser.add_argument_group('Output')
    output.add_argument('-o', dest='outdir', required=False,
                        help='Define directory where detailed output should be places - If not provided, summary will be printed to std-out')
    output.add_argument('-n', dest='outname', required=False,
                        help='Define output file name - Mandatory is -o is provided: <outname>_<contig_id>_ORF_Comparison.csv')

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-v', dest='verbose', default='False', type=eval, choices=[True, False],
                        help='Default - False: Print out runtime status')
    options = parser.parse_args()
    comparator(options)

if __name__ == "__main__":
    main()
    print("Complete")

