import collections
import sys
try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def GFF(**kwargs):
    tool_pred, genome,types = list(kwargs.values())
    GFF_ORFs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)
    with open(tool_pred, 'r') as gff_input:
        for line in gff_input:
            if '#' not in line:
                line = line.split('\t')
                gene_types = types.split(',')
                if any(gene_type == line[2] for gene_type in gene_types)and len(line) == 9:  # line[2] for normalrun
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]
                    name = line[8].split('Name=')[1].split(';')[0] # Issue with multiple records for each gene.
                    if '-' in strand:  # Reverse Compliment starts and stops adjusted
                        r_start = genome_size - stop
                        r_stop = genome_size - start
                        startCodon = genome_rev[r_start:r_start + 3]
                        stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                    elif '+' in strand:
                        startCodon = genome[start - 1:start + 2]
                        stopCodon = genome[stop - 3:stop]
                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon, line[2],name] # This needs to detect the type
                    GFF_ORFs.update({po: orf})
                # elif "CDS" in line[2]:
                #     sys.exit("SAS")

    GFF_ORFs = sortORFs(GFF_ORFs)
    return GFF_ORFs
