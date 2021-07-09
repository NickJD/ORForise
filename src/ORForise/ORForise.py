import argparse
import collections
import csv
import sys
from importlib import import_module
from Comparator import tool_comparison


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tool', required=True, help='Which tool to analyse?')
parser.add_argument('-tp', '--tool_pred', required=True, help='Tool genome prediction file - Different Tool Parameters'
                                                               ' are compared induvidually via separate files')
parser.add_argument('-dna', '--genome_DNA', required=True, help='Genome DNA file')
parser.add_argument('-gff', '--genome_GFF', required=True, help='Genome GFF file')
parser.add_argument('-o', '--outname', required=False, help='Define Output filename -- Default: Tool_Name/Genome_Name.csv')
args = parser.parse_args()

def comparator(tool,tool_pred,genome_DNA,genome_GFF,outname):
    genome_Seq = ""
    with open(genome_DNA, 'r') as genome:
        for line in genome:
            line = line.replace("\n","")
            if not line.startswith('>'):
                genome_Seq += str(line)
    ##############################################
    genes = collections.OrderedDict() # Order is important
    count = 0
    with open(genome_GFF,'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "CDS" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]
                    gene = str(start) + ',' + str(stop) + ',' + strand
                    genes.update({count: gene})
                    count +=1
            except IndexError:
                continue
    #############################################
    try:
        tool_predictions = import_module('ORForise.Tools.'+tool+'.'+tool, package='my_current_pkg')
        tool_predictions = getattr(tool_predictions,tool)
    except ModuleNotFoundError:
        sys.exit("Tool not available")
    orfs = tool_predictions(tool_pred,genome_Seq)
    all_Metrics, all_rep_Metrics, start_precision, stop_precision,other_starts, other_stops, perfect_Matches, missed_genes, unmatched_orfs, undetected_gene_metrics, unmatched_orf_metrics, gene_coverage_genome, multi_Matched_ORFs, partial_Hits = tool_comparison(genes,orfs,genome_Seq)
    ############################################# To get default output filename from input file details
    genome_name = genome_GFF.split('/')[-1].split('.')[0]
    if not outname:
        out = tool + '_' + genome_name
        outname = "Tools/" + tool + '/' + out + '.csv'
    metric_description = list(all_Metrics.keys())
    metrics = list(all_Metrics.values())
    rep_metric_description = list(all_rep_Metrics.keys())
    rep_metrics = list(all_rep_Metrics.values())
    with open(outname, 'w', newline='\n', encoding='utf-8') as out_file: # Clear write out of report
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
        tool_out.writerow(['ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'])
        tool_out.writerow(undetected_gene_metrics)
        ####
        tool_out.writerow(['Perfect_Match_Genes:'])
        for key,value in perfect_Matches.items():
            key = key.split(',')
            id = ('>' + genome_name + '_' + key[0] + '_' + key[1] + '_' + key[2])
            tool_out.writerow([id + '\n' + value + '\n'])
        ####
        tool_out.writerow(['\nUndetected_Genes:'])
        for key,value in missed_genes.items():
            key = key.split(',')
            id = ('>' + genome_name + '_' + key[0] + '_' + key[1] + '_' + key[2])
            tool_out.writerow([id + '\n' + value + '\n'])
        tool_out.writerow(['\nORFs_Without_Corresponding_Gene_In_Ensembl_Metrics:'])
        tool_out.writerow(['ATG_Start,GTG_Start,TTG_Start,ATT_Start,CTG_Start,Alternative_Start_Codon,TGA_Stop,TAA_Stop,TAG_Stop,Alternative_Stop_Codon,Median_Length,ORFs_on_Positive_Strand,ORFs_on_Negative_Strand'])
        tool_out.writerow(unmatched_orf_metrics)
        tool_out.writerow(['ORF_Without_Corresponding_Gene_in_Ensembl:'])
        for key, value in unmatched_orfs.items():
            key = key.split(',')
            id = ('>'+tool+'_'+key[0]+'_'+key[1]+'_'+key[2])
            tool_out.writerow([id + '\n' + value])
        tool_out.writerow(['\nORFs_Which_Detected_more_than_one_Gene:'])

        try:
            for key, value in multi_Matched_ORFs.items():
                key = key.split(',')
                multi = ('ORF:'+key[0]+'-'+key[1]+'_Genes:'+'|'.join(value))
                tool_out.writerow([multi])
        except IndexError:
            pass

        tool_out.writerow(['\n\nPartial_Gene_Hits:'])
        for key, seqs in partial_Hits.items():
            key = key.split(';')
            gene_Seq = seqs[0]
            orf_Seq = seqs[1]
            partial = (key[0]+'\n'+gene_Seq+'\n'+key[1]+'\n'+orf_Seq+'\n')
            tool_out.writerow([partial])




if __name__ == "__main__":
    comparator(**vars(args))

    print("Complete")
