from importlib import import_module

import argparse
import collections
import csv

from Comparator import tool_comparison

###################


def comparator(tool, input_to_analyse, storfs_to_find_missing, genome_to_compare):
    genome_Seq = ""
    with open('Genomes/' + genome_to_compare + '.fa', 'r') as genome:
        for line in genome:
            line = line.replace("\n", "")
            if ">" not in line:
                genome_Seq += str(line)
    ##############################################
    genes = collections.OrderedDict()
    count = 0
    with open('Tools/StORF_Undetected/' + input_to_analyse, 'r') as genome_gff:  # Get list of missed genes
        for line in genome_gff:
            if ">" in line:
                line = line.strip()
                Start = int(line.split('_')[1])
                Stop = int(line.split('_')[2])
                Strand = line.split('_')[3]
                Gene = str(Start) + ',' + str(Stop) + ',' + Strand
                genes.update({count: Gene})
                count += 1
    ##################################
    tool_predictions = import_module('Tools.' + tool + '.' + tool)
    tool_predictions = getattr(tool_predictions, tool)
    orfs = tool_predictions(storfs_to_find_missing, genome_Seq)
    all_Metrics, all_rep_Metrics, start_precision, stop_precision, other_starts, other_stops, missed_genes, unmatched_orfs, undetected_gene_metrics, unmatched_orf_metrics, gene_coverage_genome, multi_Matched_ORFs, partial_Hits = tool_comparison(
        genes, orfs, genome_Seq)
    outname = tool + '_' + genome_to_compare
    metric_description = list(all_Metrics.keys())
    metrics = list(all_Metrics.values())
    rep_metric_description = list(all_rep_Metrics.keys())
    rep_metrics = list(all_rep_Metrics.values())
    with open("Tools/" + tool + '/' + outname + '.csv', 'w', newline='\n',
              encoding='utf-8') as out_file:  # Clear write out of report
        tool_out = csv.writer(out_file, quoting=csv.QUOTE_NONE, escapechar=" ")
        tool_out.writerow(['Representative_Metrics:'])
        tool_out.writerow(rep_metric_description)
        tool_out.writerow(rep_metrics)
        tool_out.writerow(['All_Metrics:'])
        tool_out.writerow(metric_description)
        tool_out.writerow(metrics)
        tool_out.writerow(['CDS_Gene_Coverage_of_Genome:'])
        tool_out.writerow([gene_coverage_genome])
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
        tool_out.writerow(['Undetected_Genes:'])
        for key, value in missed_genes.items():
            key = key.split(',')
            id = ('>' + genome_to_compare + '_' + key[0] + '_' + key[1] + '_' + key[2])
            tool_out.writerow([id + '\n' + value])
        tool_out.writerow(['\nORFs_Without_Corresponding_Gene_In_Ensembl_Metrics:'])
        tool_out.writerow([
            'ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'])
        tool_out.writerow(unmatched_orf_metrics)
        tool_out.writerow(['ORF_Without_Corresponding_Gene_in_Ensembl:'])
        for key, value in unmatched_orfs.items():
            key = key.split(',')
            id = ('>' + tool + '_' + key[0] + '_' + key[1] + '_' + key[2])
            tool_out.writerow([id + '\n' + value])
        tool_out.writerow(['\nORFs_Which_Detected_more_than_one_Gene:'])

        try:
            for key, value in multi_Matched_ORFs.items():
                key = key.split(',')
                value = value[1].split(',')
                multi = ('ORF:' + key[0] + '-' + key[1] + '_Gene:' + value[0] + '-' + value[1])
                tool_out.writerow([multi])
        except IndexError:
            pass

        tool_out.writerow(['\n\nPartial_Gene_Hits:'])
        for key, seqs in partial_Hits.items():
            key = key.split(';')
            gene_Seq = seqs[0]
            orf_Seq = seqs[1]
            partial = (key[0] + '\n' + gene_Seq + '\n' + key[1] + '\n' + orf_Seq + '\n')
            tool_out.writerow([partial])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tool', default='GFF', help='Which tool/format would you analyse with StORF-R?')
    parser.add_argument('-i', '--input_to_analyse', default='', help='Location of file containing missed genes')
    parser.add_argument('-stf', '--storfs_to_find_missing', default='', help='STORFs to find missing.')
    parser.add_argument('-g', '--genome_to_compare', default='', help='Which genome to analyse?')
    args = parser.parse_args()

    comparator(**vars(args))

if __name__ == "__main__":
    main()
