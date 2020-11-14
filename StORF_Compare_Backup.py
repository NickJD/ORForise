import collections
import csv
import argparse
from importlib import import_module
from Comparator import tool_comparison

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tool', default='GFF', help='Which tool/format would you like to analyse?')
parser.add_argument('-i', '--input_to_analyse', default='', help='Location of file to compare.')
parser.add_argument('-stf', '--storfs_to_find_missing', default='', help='STORFs to find missing.')
parser.add_argument('-g', '--genome_to_compare', default='', help='Which genome to analyse?')
args = parser.parse_args()

def comparator(tool,input_to_analyse,storfs_to_find_missing,genome_to_compare):
    genome_Seq = ""
    with open('genomes/'+genome_to_compare+'.fa', 'r') as genome:
        for line in genome:
            line = line.replace("\n","")
            if ">" not in line:
                genome_Seq += str(line)
    ##############################################
    genes = collections.OrderedDict()
    count = 0
    with open('Tools/StORF_Undetected/'+input_to_analyse,'r') as genome_gff:
        for line in genome_gff:
            if ">" in line:
                line = line.strip()
                Start = int(line.split('_')[1])
                Stop = int(line.split('_')[2])
                Strand = line.split('_')[3]
                Gene = str(Start) + ',' + str(Stop) + ',' + Strand
                genes.update({count: Gene})
                count +=1
# ##################################
    tool_predictions = import_module('Tools.'+tool+'.'+tool)
    tool_predictions = getattr(tool_predictions,tool)
    orfs = tool_predictions(storfs_to_find_missing,genome_Seq)
    metric_description,metrics, rep_metric_description, rep_metrics, start_precision, stop_precision,other_starts, other_stops, missed_genes, unmatched_orfs, undetected_gene_metrics, unmatched_orf_metrics, gene_coverage_genome = tool_comparison(genes,orfs,genome_Seq)

    outname = input_to_analyse.split('.')[0]

    with open("Tools/" + tool + '/' + outname + '.csv', 'w', newline='\n',
              encoding='utf-8') as out_file:  # Clear write out of report
        tool_out = csv.writer(out_file, quoting=csv.QUOTE_NONE, escapechar=" ")
        tool_out.writerow(['Representative Metrics:'])
        tool_out.writerow(rep_metric_description)
        tool_out.writerow(rep_metrics)
        tool_out.writerow(['All Metrics:'])
        tool_out.writerow(metric_description)
        tool_out.writerow(metrics)
        tool_out.writerow(['CDS Gene Coverage of Genome: '])
        tool_out.writerow([format(gene_coverage_genome, '.2f')])
        tool_out.writerow(['Start Position Difference:'])
        tool_out.writerow(start_precision)
        tool_out.writerow(['Stop Position Difference:'])
        tool_out.writerow(stop_precision)
        tool_out.writerow(['Alternative Starts Predicted:'])
        tool_out.writerow(other_starts)
        tool_out.writerow(['Alternative Stops Predicted:'])
        tool_out.writerow(other_stops)
        tool_out.writerow(['Undetected Gene Metrics:'])
        tool_out.writerow([
                              'ATG Start,GTG Start,TTG Start,ATT Start,CTG Start,Alternative Start Codon,TGA Stop,TAA Stop,TAG Stop,Alternative Stop Codon,Median Length,Genes on Positive Strand,Genes on Negative Strand'])
        tool_out.writerow(undetected_gene_metrics)
        tool_out.writerow(['Undetected Genes:'])
        for key, value in missed_genes.items():
            key = key.split(',')
            id = ('>' + genome_to_compare + '_' + key[0] + '_' + key[1] + '_' + key[2])
            tool_out.writerow([id + '\n' + value])
        tool_out.writerow(['\n\n\nORFs without corresponding gene in Ensembl Metrics:'])
        tool_out.writerow([
                              'ATG Start,GTG Start,TTG Start,ATT Start,CTG Start,Alternative Start Codon,TGA Stop,TAA Stop,TAG Stop,Alternative Stop Codon,Median Length,ORFs on Positive Strand,ORFs on Negative Strand'])
        tool_out.writerow(unmatched_orf_metrics)
        tool_out.writerow(['ORFs without corresponding gene in Ensembl:'])
        for key, value in unmatched_orfs.items():
            key = key.split(',')
            id = ('>' + tool + '_' + key[0] + '_' + key[1] + '_' + key[2])
            tool_out.writerow([id + '\n' + value])

if __name__ == "__main__":
    comparator(**vars(args))
