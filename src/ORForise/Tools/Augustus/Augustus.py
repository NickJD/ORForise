import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def Augustus(tool_pred, genome):
    augustus_ORFs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)
    with open(tool_pred, 'r') as Augustus_input:
        for line in Augustus_input:
            line = line.split()
            if len(line) == 12 and "CDS" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
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
                augustus_ORFs.update({po: orf})

    augustus_ORFs = sortORFs(augustus_ORFs)
    return augustus_ORFs
