import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def StORF_Undetected(tool_pred, genome):
    storf_orfs = collections.OrderedDict()
    genome_size = len(genome)
    genome_Rev = revCompIterative(genome)
    with open(tool_pred, 'r') as storf_input:
        for line in storf_input:
            line = line.split()
            if "StORF" in line[1] and "ORF" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                if '-' in strand:  # Reverse Compliment starts and stops adjusted
                    r_start = genome_size - stop
                    r_stop = genome_size - start
                    startCodon = genome_Rev[r_start:r_start + 3]
                    stopCodon = genome_Rev[r_stop - 2:r_stop + 1]
                elif '+' in strand:
                    startCodon = genome[start - 1:start - 1 + 3]
                    stopCodon = genome[stop - 3:stop - 1 + 1]
                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                storf_orfs.update({po: orf})

    storf_orfs = sortORFs(storf_orfs)
    return storf_orfs
