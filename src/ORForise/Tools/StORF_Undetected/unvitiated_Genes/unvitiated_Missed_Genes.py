import argparse
import collections

parser = argparse.ArgumentParser()
parser.add_argument('-mg', '--missing_genes', default='', help='Which set of genes to check?')
parser.add_argument('-g', '--genome_to_compare', required=True,
                    help='Which genome to analyse? Genome files have same prefix'
                         ' - .fa and .gff appended')
args = parser.parse_args()


def comparator(missing_genes, genome_to_compare):
    missed_genes = collections.OrderedDict()
    non_vitiated_genes = []
    with open('../' + missing_genes, 'r') as m_genes:
        for line in m_genes:
            if line.startswith('>'):
                line = line.split('_')
                start = int(line[1])
                stop = int(line[2])
                m_gene_Set = set(range(start, stop + 1))
                missed_genes.update({str(start) + '_' + str(stop): m_gene_Set})
                non_vitiated_genes.append(str(start) + '_' + str(stop))
    ##############################################

    with open('../../Prodigal/Prodigal_' + genome_to_compare + '.gff', 'r') as prodigal_input:
        for line in prodigal_input:
            line = line.split()
            if "Prodigal" in line[1] and "CDS" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                pred_set = set(range(start, stop + 1))
                for missed, g_set in missed_genes.items():
                    overlap = len(pred_set.intersection(g_set))
                    if overlap > 50:
                        print(start)
                        try:
                            non_vitiated_genes.remove(missed)
                            print(missed)
                        except ValueError:
                            continue
    print(len(non_vitiated_genes))


if __name__ == "__main__":
    comparator(**vars(args))
