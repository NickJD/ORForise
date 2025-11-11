#!/usr/bin/env python3
import collections

# Constants
SHORT_ORF_LENGTH = 300
MIN_COVERAGE = 75
ORForise_Version = 'v1.5.1'


def revCompIterative(watson):  # Gets Reverse Complement
    return watson.upper()[::-1].translate(str.maketrans("ATCGRYKMVBHD","TAGCYRMKBVDH"))


def sortORFs(tool_ORFs):  # Will only sort by given start position
    tool_ORFs_Sorted = sorted(tool_ORFs.items(), key=lambda v: int(v[0].split(",")[0]))
    tool_ORFs_Sorted = collections.OrderedDict(tool_ORFs_Sorted)
    return tool_ORFs_Sorted


def sortGenes(Genes):  # Will sort by given start position and then rearrange for given stop
    Genes_Sorted_list = sorted(Genes.values(), key=lambda v: int(v[0]))
    Genes_Sorted = []
    for idx,gene in enumerate(Genes_Sorted_list):
        Genes_Sorted.append([idx,gene])
    Genes_Sorted = collections.OrderedDict(Genes_Sorted)
    prev_stop = 0
    for pos, detail in Genes_Sorted.items():
        if detail[1] < prev_stop:
            Genes_Sorted[pos], Genes_Sorted[pos-1] = Genes_Sorted[pos-1], Genes_Sorted[pos]
        prev_stop = detail[1]
    return Genes_Sorted


def gff_load(options,gff_in,dna_regions):
    count = 0
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split('\t')
        if line.startswith('\n') or line.startswith('#') or 'European Nucleotide Archive' in line:  # Not to crash on empty lines in GFF
            continue
        elif options.gene_ident[0] == 'ID=gene':
            if line_data[0] in dna_regions and options.gene_ident[0] in line_data[8]:
                start = int(line_data[3])
                stop = int(line_data[4])
                strand = line_data[6]
                gene_details = [start,stop,strand]
                dna_regions[line_data[0]][2].append({count:gene_details}) # This will add to list
                count += 1
        else:
            try:
                if line_data[2] == 'region':
                    continue
                elif line_data[0] in dna_regions:
                    if any(gene_type in line_data[2] for gene_type in options.gene_ident): # line[2] for normal run
                        start = int(line_data[3])
                        stop = int(line_data[4])
                        strand = line_data[6]
                        gene_details = [start, stop, strand]
                        if gene_details not in dna_regions[line_data[0]][2]:
                            dna_regions[line_data[0]][2].append({count:gene_details}) # This will add to list
                            count += 1
            except IndexError:
                continue
    return dna_regions


def fasta_load(fasta_in):
    dna_regions = collections.OrderedDict()
    first = True
    if '>' in fasta_in.readline().rstrip():
        fasta_in.seek(0)
        #### Default for when presented with standard fasta file
        for line in fasta_in:
            line = line.strip()
            if line.startswith('>') and first == False:  # Check if first seq in file
                dna_region_length = len(seq)
                dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
                seq = ''
                dna_region_id = line.split()[0].replace('>', '')
            elif line.startswith('>'):
                seq = ''
                dna_region_id = line.split()[0].replace('>', '')
            else:
                seq += str(line)
                first = False
        dna_region_length = len(seq)
        dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
    elif '##' in  fasta_in.readline().rstrip(): # Clunky and may fall over
        fasta_in.seek(0)
        #### Called when presented with Prokka GFF file so must get fasta from inside it
        ### Get to genome seq
        at_FASTA = False
        for line in fasta_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
            if line.startswith('##FASTA'):  # Not to crash on empty lines in GFF
                at_FASTA = True
            elif at_FASTA == True:
                line = line.strip()
                if line.startswith('>') and first == False:  # Check if first seq in file
                    dna_region_length = len(seq)
                    dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
                    seq = ''
                    dna_region_id = line.split()[0].replace('>', '')
                elif line.startswith('>'):
                    seq = ''
                    dna_region_id = line.split()[0].replace('>', '')
                else:
                    seq += str(line)
                    first = False
        dna_region_length = len(seq)
        dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})

    return dna_regions


def get_rep_metrics(result):
    rep_metric_description = ('Percentage_of_Genes_Detected,Percentage_of_ORFs_that_Detected_a_Gene,'
                              'Percent_Difference_of_All_ORFs,Median_Length_Difference,Percentage_of_Perfect_Matches,'
                              'Median_Start_Difference_of_Matched_ORFs,Median_Stop_Difference_of_Matched_ORFs,'
                              'Percentage_Difference_of_Matched_Overlapping_CDSs,Percent_Difference_of_Short-Matched-ORFs,'
                              'Precision,Recall,False_Discovery_Rate')
    rep_metrics = [result['rep_metrics']['Percentage_of_Genes_Detected'],
                   result['rep_metrics']['Percentage_of_ORFs_that_Detected_a_Gene'],
                   result['rep_metrics']['Percent_Difference_of_All_ORFs'],
                   result['rep_metrics']['Median_Length_Difference'],
                   result['rep_metrics']['Percentage_of_Perfect_Matches'],
                   result['rep_metrics']['Median_Start_Difference_of_Matched_ORFs'],
                   result['rep_metrics']['Median_Stop_Difference_of_Matched_ORFs'],
                   result['rep_metrics']['Percentage_Difference_of_Matched_Overlapping_CDSs'],
                   result['rep_metrics']['Percent_Difference_of_Short-Matched-ORFs'],
                   result['rep_metrics']['Precision'],
                   result['rep_metrics']['Recall'],
                   result['rep_metrics']['False_Discovery_Rate']]
    return rep_metric_description, rep_metrics


def get_all_metrics(result):
    all_metric_description = ('Number_of_ORFs,Percent_Difference_of_All_ORFs,Number_of_ORFs_that_Detected_a_Gene,'
         'Percentage_of_ORFs_that_Detected_a_Gene,Number_of_Genes_Detected,Percentage_of_Genes_Detected,'
         'Median_Length_of_All_ORFs,Median_Length_Difference,Minimum_Length_of_All_ORFs,Minimum_Length_Difference,'
         'Maximum_Length_of_All_ORFs,Maximum_Length_Difference,Median_GC_content_of_All_ORFs,'
         'Percent_Difference_of_All_ORFs_Median_GC,Median_GC_content_of_Matched_ORFs,'
         'Percent_Difference_of_Matched_ORF_GC,Number_of_ORFs_which_Overlap_Another_ORF,'
         'Percent_Difference_of_Overlapping_ORFs,Maximum_ORF_Overlap,Median_ORF_Overlap,'
         'Number_of_Matched_ORFs_Overlapping_Another_ORF,Percentage_Difference_of_Matched_Overlapping_CDSs,'
         'Maximum_Matched_ORF_Overlap,Median_Matched_ORF_Overlap,Number_of_Short-ORFs,Percent_Difference_of_Short-ORFs,'
         'Number_of_Short-Matched-ORFs,Percent_Difference_of_Short-Matched-ORFs,Number_of_Perfect_Matches,'
         'Percentage_of_Perfect_Matches,Number_of_Perfect_Starts,Percentage_of_Perfect_Starts,Number_of_Perfect_Stops,'
         'Percentage_of_Perfect_Stops,Number_of_Out_of_Frame_ORFs,Number_of_Matched_ORFs_Extending_a_Coding_Region,'
         'Percentage_of_Matched_ORFs_Extending_a_Coding_Region,Number_of_Matched_ORFs_Extending_Start_Region,'
         'Percentage_of_Matched_ORFs_Extending_Start_Region,Number_of_Matched_ORFs_Extending_Stop_Region,'
         'Percentage_of_Matched_ORFs_Extending_Stop_Region,Number_of_All_ORFs_on_Positive_Strand,'
         'Percentage_of_All_ORFs_on_Positive_Strand,Number_of_All_ORFs_on_Negative_Strand,'
         'Percentage_of_All_ORFs_on_Negative_Strand,Median_Start_Difference_of_Matched_ORFs,'
         'Median_Stop_Difference_of_Matched_ORFs,ATG_Start_Percentage,GTG_Start_Percentage,TTG_Start_Percentage,'
         'ATT_Start_Percentage,CTG_Start_Percentage,Other_Start_Codon_Percentage,TAG_Stop_Percentage,'
         'TAA_Stop_Percentage,TGA_Stop_Percentage,Other_Stop_Codon_Percentage,True_Positive,False_Positive,'
         'False_Negative,Precision,Recall,False_Discovery_Rate,Nucleotide_True_Positive,Nucleotide_False_Positive,'
         'Nucleotide_True_Negative,Nucleotide_False_Negative,Nucleotide_Precision,Nucleotide_Recall,'
         'Nucleotide_False_Discovery_Rate,ORF_Nucleotide_Coverage_of_Genome,Matched_ORF_Nucleotide_Coverage_of_Genome')
    all_metrics = rep_metrics = [result['pred_metrics']['Number_of_ORFs'],
                                    result['pred_metrics']['Percent_Difference_of_All_ORFs'],
                                    result['pred_metrics']['Number_of_ORFs_that_Detected_a_Gene'],
                                    result['pred_metrics']['Percentage_of_ORFs_that_Detected_a_Gene'],
                                    result['pred_metrics']['Number_of_Genes_Detected'],
                                    result['pred_metrics']['Percentage_of_Genes_Detected'],
                                    result['pred_metrics']['Median_Length_of_All_ORFs'],
                                    result['pred_metrics']['Median_Length_Difference'],
                                    result['pred_metrics']['Minimum_Length_of_All_ORFs'],
                                    result['pred_metrics']['Minimum_Length_Difference'],
                                    result['pred_metrics']['Maximum_Length_of_All_ORFs'],
                                    result['pred_metrics']['Maximum_Length_Difference'],
                                    result['pred_metrics']['Median_GC_content_of_All_ORFs'],
                                    result['pred_metrics']['Percent_Difference_of_All_ORFs_Median_GC'],
                                    result['pred_metrics']['Median_GC_content_of_Matched_ORFs'],
                                    result['pred_metrics']['Percent_Difference_of_Matched_ORF_GC'],
                                    result['pred_metrics']['Number_of_ORFs_which_Overlap_Another_ORF'],
                                    result['pred_metrics']['Percent_Difference_of_Overlapping_ORFs'],
                                    result['pred_metrics']['Maximum_ORF_Overlap'],
                                    result['pred_metrics']['Median_ORF_Overlap'],
                                    result['pred_metrics']['Number_of_Matched_ORFs_Overlapping_Another_ORF'],
                                    result['pred_metrics']['Percentage_Difference_of_Matched_Overlapping_CDSs'],
                                    result['pred_metrics']['Maximum_Matched_ORF_Overlap'],
                                    result['pred_metrics']['Median_Matched_ORF_Overlap'],
                                    result['pred_metrics']['Number_of_Short-ORFs'],
                                    result['pred_metrics']['Percent_Difference_of_Short-ORFs'],
                                    result['pred_metrics']['Number_of_Short-Matched-ORFs'],
                                    result['pred_metrics']['Percent_Difference_of_Short-Matched-ORFs'],
                                    result['pred_metrics']['Number_of_Perfect_Matches'],
                                    result['pred_metrics']['Percentage_of_Perfect_Matches'],
                                    result['pred_metrics']['Number_of_Perfect_Starts'],
                                    result['pred_metrics']['Percentage_of_Perfect_Starts'],
                                    result['pred_metrics']['Number_of_Perfect_Stops'],
                                    result['pred_metrics']['Percentage_of_Perfect_Stops'],
                                    result['pred_metrics']['Number_of_Out_of_Frame_ORFs'],
                                    result['pred_metrics']['Number_of_Matched_ORFs_Extending_a_Coding_Region'],
                                    result['pred_metrics']['Percentage_of_Matched_ORFs_Extending_a_Coding_Region'],
                                    result['pred_metrics']['Number_of_Matched_ORFs_Extending_Start_Region'],
                                    result['pred_metrics']['Percentage_of_Matched_ORFs_Extending_Start_Region'],
                                    result['pred_metrics']['Number_of_Matched_ORFs_Extending_Stop_Region'],
                                    result['pred_metrics']['Percentage_of_Matched_ORFs_Extending_Stop_Region'],
                                    result['pred_metrics']['Number_of_All_ORFs_on_Positive_Strand'],
                                    result['pred_metrics']['Percentage_of_All_ORFs_on_Positive_Strand'],
                                    result['pred_metrics']['Number_of_All_ORFs_on_Negative_Strand'],
                                    result['pred_metrics']['Percentage_of_All_ORFs_on_Negative_Strand'],
                                    result['pred_metrics']['Median_Start_Difference_of_Matched_ORFs'],
                                    result['pred_metrics']['Median_Stop_Difference_of_Matched_ORFs'],
                                    result['pred_metrics']['ATG_Start_Percentage'],
                                    result['pred_metrics']['GTG_Start_Percentage'],
                                    result['pred_metrics']['TTG_Start_Percentage'],
                                    result['pred_metrics']['ATT_Start_Percentage'],
                                    result['pred_metrics']['CTG_Start_Percentage'],
                                    result['pred_metrics']['Other_Start_Codon_Percentage'],
                                    result['pred_metrics']['TAG_Stop_Percentage'],
                                    result['pred_metrics']['TAA_Stop_Percentage'],
                                    result['pred_metrics']['TGA_Stop_Percentage'],
                                    result['pred_metrics']['Other_Stop_Codon_Percentage'],
                                    result['pred_metrics']['True_Positive'],
                                    result['pred_metrics']['False_Positive'],
                                    result['pred_metrics']['False_Negative'],
                                    result['pred_metrics']['Precision'],
                                    result['pred_metrics']['Recall'],
                                    result['pred_metrics']['False_Discovery_Rate'],
                                    result['pred_metrics']['Nucleotide_True_Positive'],
                                    result['pred_metrics']['Nucleotide_False_Positive'],
                                    result['pred_metrics']['Nucleotide_True_Negative'],
                                    result['pred_metrics']['Nucleotide_False_Negative'],
                                    result['pred_metrics']['Nucleotide_Precision'],
                                    result['pred_metrics']['Nucleotide_Recall'],
                                    result['pred_metrics']['Nucleotide_False_Discovery_Rate'],
                                    result['pred_metrics']['ORF_Nucleotide_Coverage_of_Genome'],
                                    result['pred_metrics']['Matched_ORF_Nucleotide_Coverage_of_Genome']]


    return all_metric_description, all_metrics