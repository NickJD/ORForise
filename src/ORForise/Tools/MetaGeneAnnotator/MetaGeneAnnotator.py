import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def MetaGeneAnnotator(tool_pred, genome):
    metaGeneAnnotator_ORFs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)
    with open(tool_pred, 'r') as MetaGeneAnnotator_input:
        for line in MetaGeneAnnotator_input:
            line = line.split()
            if len(line) == 11:
                if "gene_" in line[0]:
                    start = int(line[1])
                    stop = int(line[2])
                    strand = line[3]
                    if '-' in strand:  # Reverse Compliment starts and stops adjusted
                        r_start = genome_size - stop
                        r_stop = genome_size - start
                        startCodon = genome_rev[r_start:r_start + 3]
                        stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                    elif '+' in strand:
                        startCodon = genome[start - 1:start + 2]
                        stopCodon = genome[stop - 3:stop]
                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon, 'CDS']
                    metaGeneAnnotator_ORFs.update({po: orf})

    metaGeneAnnotator_ORFs = sortORFs(metaGeneAnnotator_ORFs)
    return metaGeneAnnotator_ORFs
