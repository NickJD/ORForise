#!/usr/bin/env python3
import collections

# Constants
SHORT_ORF_LENGTH = 300
MIN_COVERAGE = 75
ORForise_Version = 'v1.4.2'


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