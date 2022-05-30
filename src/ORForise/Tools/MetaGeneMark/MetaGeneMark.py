import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def MetaGeneMark(tool_pred, genome):
    metaGeneMarkORFs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)
    with open(tool_pred, 'r') as metaGeneMark_input:
        for line in metaGeneMark_input:
            line = line.split()
            if len(line) == 19:
                if 'GeneMark.hmm' in line[4] and "CDS" in line[5]:
                    start = int(line[6])
                    stop = int(line[7])
                    strand = line[9]
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
                    metaGeneMarkORFs.update({po: orf})

    metaGeneMarkORFs = sortORFs(metaGeneMarkORFs)
    return metaGeneMarkORFs
