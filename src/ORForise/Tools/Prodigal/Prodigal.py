import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def Prodigal(*args):
    tool_pred = args[0]
    genome = args[1]
    #types = args[2]
    prodigalORFs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)
    with open(tool_pred, 'r') as prodigal_input:
        for line in prodigal_input:
            line = line.split()
            if "Prodigal" in line[1] and "CDS" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                info = line[8]
                if '-' in strand:  # Reverse Compliment starts and stops adjusted
                    r_start = genome_size - stop
                    r_stop = genome_size - start
                    startCodon = genome_rev[r_start:r_start + 3]
                    stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                elif '+' in strand:
                    startCodon = genome[start - 1:start + 2]
                    stopCodon = genome[stop - 3:stop]
                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon, 'CDS', 'Prodigal|'+info]
                prodigalORFs.update({po: orf})

    prodigalORFs = sortORFs(prodigalORFs)
    return prodigalORFs
