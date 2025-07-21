import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def GLIMMER_3(*args):
    tool_pred = args[0]
    dna_regions = args[1]
    GLIMMER_ORFs = collections.OrderedDict()
    for dna_region in dna_regions:
        GLIMMER_ORFs[dna_region] = collections.OrderedDict()
    for dna_region in dna_regions:
        genome = dna_regions[dna_region][0]
        genome_size = len(genome)
        genome_rev = revCompIterative(genome)
        with open(tool_pred,
                  'r') as glimmer_input:  # GLIMMER_3 reverses the start and stop positions for ORFS on the negative strand
            for line in glimmer_input:
                if '>' not in line:  # This will not work with multiple contigs
                    line = line.split()
                    if len(line) == 5 and "orf" in line[0] and dna_region in line[0]:
                        if '-' in line[3]:  # Reverse Compliment starts and stops adjusted -  Switched to match Sense Strand
                            start = int(line[2])
                            stop = int(line[1])
                            strand = '-'
                            r_start = genome_size - stop
                            r_stop = genome_size - start
                            startCodon = genome_rev[r_start:r_start + 3]
                            stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                        elif '+' in line[3]:
                            start = int(line[1])
                            stop = int(line[2])
                            strand = '+'
                            startCodon = genome[start - 1:start + 3]
                            stopCodon = genome[stop - 3:stop]
                        po = str(start) + ',' + str(stop)
                        orf = [strand, startCodon, stopCodon, 'CDS', 'GLIMMER_3']
                        GLIMMER_ORFs.update({po: orf})

    for group in GLIMMER_ORFs:
        GLIMMER_ORFs[group] = sortORFs(GLIMMER_ORFs[group])
    return GLIMMER_ORFs
