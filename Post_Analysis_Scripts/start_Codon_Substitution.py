import argparse
from ORForise.Tools.utils import *
import collections


parser = argparse.ArgumentParser()
parser.add_argument('-r', '--results', required=True, help='Which output to look at?')
parser.add_argument('-a', '--annotation', required=True, help='Genome Annotation File')

args = parser.parse_args()


def revCompIterative(watson):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev:
        crick += complements[nt]
    return crick


def partial_gene_read(results,partial_genes):
    #partial Genes Read-In
    orf_Lengths = []
    gene_Lengths = []
    results_in = open('../Tools/'+results)
    read = False
    prev = ''
    for line in results_in:
        line = line.strip()
        if read == True:
            if line.startswith('Gene:'):
                line = line.replace('Gene:','')
                entry = line.split('_')
                g_Pos = entry[0]+'_'+entry[1]
                strand = entry[2]
                prev = 'Gene'
            elif line.startswith('ORF:'):
                line = line.replace('ORF:', '')
                entry = line.split('_')
                o_Pos = entry[0] + '_' + entry[1]
                prev = 'ORF'
            elif line:
                if 'Gene' in prev:
                    g_Seq = line.strip()
                    gene_Lengths.append(len(g_Seq))
                elif 'ORF' in prev:
                    o_Seq = line.strip()
                    orf_Lengths.append(len(o_Seq))
            elif not line:
                partial_genes.update({g_Pos:[strand,g_Seq,o_Pos,o_Seq]})
        if line.startswith('Partial_Gene_Hits:'):
            read = True

    return partial_genes,orf_Lengths,gene_Lengths


def detail_transfer(genes,partial_genes):
    for partial,m_details in partial_genes.items():
        try:
            details = genes[partial]
            gc = details[2]
            up_Overlap = details[3]
            down_Overlap = details[4]
            m_details.insert(2,gc)
            m_details.insert(3,up_Overlap)
            m_details.insert(4,down_Overlap)
        except KeyError:
            pass
    return partial_genes


def result_compare(results,annotation):
    genome = ""
    with open('../Genomes/' + annotation + '.fa', 'r') as genome_file:
        for line in genome_file:
            line = line.replace("\n", "")
            if ">" not in line:
                genome += str(line)

    partial_genes = collections.OrderedDict()
    partial_genes,orf_Lengths,gene_Lengths = partial_gene_read(results,partial_genes)

    strands = collections.defaultdict(int, {'-': 0, '+': 0})

    start_Codon_Substitution = collections.OrderedDict({'ATG-ATG':0,'CTG-ATG':0,'GTG-ATG':0,'TTG-ATG':0,'Alt-ATG':0,
                                                        'ATG-CTG':0,'CTG-CTG':0,'GTG-CTG':0,'TTG-CTG':0,'Alt-CTG':0,
                                                        'ATG-GTG':0,'CTG-GTG':0,'GTG-GTG':0,'TTG-GTG':0,'Alt-GTG':0,
                                                        'ATG-TTG':0,'CTG-TTG':0,'GTG-TTG':0,'TTG-TTG':0,'Alt-TTG':0,
                                                        'ATG-Alt':0,'CTG-Alt':0,'GTG-Alt':0,'TTG-Alt':0,'Alt-Alt':0})

    codon_set = ['ATG','CTG','GTG','TTG']
    for gene, data in partial_genes.items():
        strands[data[0]] +=1
        gene_Start = [data[1][0:3]]
        orf_Start = [data[3][0:3]]
        if gene_Start[0] in codon_set:
            gene_Start = gene_Start[0]
        else:
            print('Gene_Codon_Alternative:'+str(gene_Start[0]))
            gene_Start = 'Alt'
        if orf_Start[0] in codon_set:
            orf_Start = orf_Start[0]
        else:
            print('ORF_Codon_Alternative:'+str(orf_Start[0]))
            orf_Start = 'Alt'

        matrix_index = gene_Start+'-'+orf_Start
        start_Codon_Substitution[matrix_index] +=1

    ####### HERE - Need to flip the data  - GS along the top
    subs = start_Codon_Substitution.values()
    subs = list(subs)
    for i in [subs[c:c + 5] for c in range(0, len(subs), 5) if c % 5 == 0]:
        print(*i)




if __name__ == "__main__":
    result_compare(**vars(args))



