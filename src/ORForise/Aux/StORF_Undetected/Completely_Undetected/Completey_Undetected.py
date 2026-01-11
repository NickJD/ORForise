import argparse

from orderedset import OrderedSet

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--undetected_Genes', default='', help='Undected Genes.')
parser.add_argument('-t', '--tool', default='', help='Tool Used.')
args = parser.parse_args()


def un_Genes(undetected_Genes, tool_ORFs):
    count = 0
    genes = []
    with open(undetected_Genes, 'r') as undected_GFF:
        for line in undected_GFF:
            if ">" in line:
                line = line.split('_')
                g_start = int(line[1])
                g_stop = int(line[2])
                gene_ORF = str(g_start) + ',' + str(g_stop)
                genes.append(gene_ORF)
                gene_Set = OrderedSet(range(g_start, g_stop + 1))
                for t_ORF in tool_ORFs:
                    t_ORF = t_ORF.split(',')
                    t_start = int(t_ORF[0])
                    t_stop = int(t_ORF[1])
                    tool_Set = OrderedSet(range(t_start - 0, t_stop + 1))
                    overlap = len(tool_Set.intersection(gene_Set))
                    if overlap >= 1:
                        print(overlap)
                        print(line)
                        count += 1
                        break
    print(
        "Number of Undetected Genes: " + str(len(genes)) + " Number of Spoiled Undetected Genes By Tool: " + str(count))


def tool(tool, undetected_Genes):
    tool_ORFs = []
    with open(tool, 'r') as input:
        for line in input:
            line = line.split()
            if "Prodigal" in line[1] and "CDS" in line[2]:  # Modify For Tool
                t_start = int(line[3])
                t_stop = int(line[4])
                tool_ORF = str(t_start) + ',' + str(t_stop)
                tool_ORFs.append(tool_ORF)

    un_Genes(undetected_Genes, tool_ORFs)


if __name__ == "__main__":
    tool(**vars(args))
