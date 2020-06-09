import collections
import csv
import argparse
from importlib import import_module
from Comparator import tool_comparison

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tool', required=True, help='Which tool to compare?')
parser.add_argument('-i', '--input_to_analyse', required=True, help='Location of tool output to compare.')
parser.add_argument('-g', '--genome_to_compare', required=True, help='Which genome to analyse?')
args = parser.parse_args()

def comparator(tool,input_to_analyse,genome_to_compare):

    global genome_seq
    genome_seq = ""
    with open('genomes/'+genome_to_compare+'_dna.fa', 'r') as genome:
        for line in genome:
            line = line.replace("\n","")
            if ">" not in line:
                genome_seq += str(line)
    ##############################################
    genes = collections.OrderedDict()
    count = 0
    with open('genomes/'+genome_to_compare+'.gff','r') as genome_gff: # Should work for GFF3
        for line in genome_gff:
            line = line.split('\t')
            if "CDS" in line and len(line) == 9:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                gene = str(start) + ',' + str(stop) + ',' + strand
                genes.update({count: gene})
                count +=1
    #############################################
    tool_predictions = import_module('Tools.'+tool+'.'+tool)
    tool_predictions = getattr(tool_predictions,tool)
    orfs = tool_predictions(input_to_analyse,genome_seq)
    metric_description,metrics, rep_metric_description, rep_metrics, start_precision, stop_precision,other_starts, other_stops, missed_genes, unmatched_orfs, missed_gene_metrics, unmatched_orf_metrics, gene_coverage_genome = tool_comparison(genes,orfs,genome_seq)

    outname = input_to_analyse.split('.')[0]

    with open("Tools/"+tool+'/'+outname+'.csv', 'w') as out_file: # Clear write out of report
        tool_out = csv.writer(out_file, quoting=csv.QUOTE_NONE, escapechar=" ")
        tool_out.writerow(['Representative Metric Descriptions'])
        tool_out.writerow(rep_metric_description)
        tool_out.writerow(['Representative Metrics:'])
        tool_out.writerow(rep_metrics)
        tool_out.writerow(['All Metric Descriptions:'])
        tool_out.writerow(metric_description)
        tool_out.writerow(['Metrics:'])
        tool_out.writerow(metrics)
        tool_out.writerow(['CDS Gene Coverage of Genome: '])
        tool_out.writerow([format(gene_coverage_genome,'.2f')])
        tool_out.writerow(['Start Codon Difference:'])
        tool_out.writerow(start_precision)
        tool_out.writerow(['Stop Codon Difference:'])
        tool_out.writerow(stop_precision)
        tool_out.writerow(['Alternative Starts:'])
        tool_out.writerow(other_starts)
        tool_out.writerow(['Alternative Stops:'])
        tool_out.writerow(other_stops)
        tool_out.writerow(['Mised Genes: Start-Stop-Strand-Start_Codon-Stop_Codon'])
        tool_out.writerow(['ATG Start,GTG Start,TTG Start,ATT Start,CTG Start,Alternative Start Codon,TGA Stop,TAA Stop,TAG Stop,Alternative Stop Codon,Mean Length,Genes on Positive Strand,Genes on Negative Strand'])
        tool_out.writerow(missed_gene_metrics)
        for key,value in missed_genes.items():
            key = key.split(',')
            id = ('>' + genome_to_compare + '_' + key[0] + '_' + key[1] + '_' + key[2])
            tool_out.writerow([id + '\n' + value])

        tool_out.writerow(['\n\n\nUnmatched ORFs: Start-Stop-Strand-Start_Codon-Stop_Codon'])
        tool_out.writerow(['ATG Start,GTG Start,TTG Start,ATT Start,CTG Start,Alternative Start Codon,TGA Stop,TAA Stop,TAG Stop,Alternative Stop Codon,Mean Length,ORFs on Positive Strand,ORFs on Negative Strand'])
        tool_out.writerow(unmatched_orf_metrics)
        for key, value in unmatched_orfs.items():
            key = key.split(',')
            id = ('>'+tool+'_'+key[0]+'_'+key[1]+'_'+key[2])
            tool_out.writerow([id + '\n' + value])


if __name__ == "__main__":
    comparator(**vars(args))








