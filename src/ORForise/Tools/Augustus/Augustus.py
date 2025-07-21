import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def Augustus(*args):
    tool_pred = args[0]
    dna_regions = args[1]
    augustus_ORFs = collections.OrderedDict()
    for dna_region in dna_regions:
        augustus_ORFs[dna_region] = collections.OrderedDict()
    for dna_region in dna_regions:
        genome = dna_regions[dna_region][0]
        genome_size = len(genome)
        genome_rev = revCompIterative(genome)
        with open(tool_pred, 'r') as Augustus_input:
            for line in Augustus_input:
                line = line.split()
                if len(line) == 12 and dna_region in line[0] and "CDS" in line[2]:
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
                    orf = [strand, startCodon, stopCodon, 'CDS', 'Augustus']
                    augustus_ORFs.update({po: orf})

        for group in augustus_ORFs:
            augustus_ORFs[group] = sortORFs(augustus_ORFs[group])
        return augustus_ORFs
