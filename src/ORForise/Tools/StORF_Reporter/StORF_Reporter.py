import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def StORF_Reporter(tool_pred, genome):
    storf_orfs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)
    with open(tool_pred, 'r') as storf_input:
        for line in storf_input:
            if '#' not in line:
                line = line.split()
                if 'StORF_Reporter' in line[1] or 'StoRF_Reporter' in line[1]: # need to harmonise this.
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]
                    if '-' in strand:  # Reverse Compliment starts and stops adjusted
                        r_start = genome_size - stop
                        r_stop = genome_size - start
                        startCodon = genome_rev[r_start:r_start + 3]
                        stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                    elif '+' in strand:
                        startCodon = genome[start:start + 3]
                        stopCodon = genome[stop - 3:stop]
                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon, line[2]] # StORF/Con-StORF or CDS??
                    storf_orfs.update({po: orf})

    storf_orfs = sortORFs(storf_orfs)
    return storf_orfs
