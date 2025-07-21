import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def FGENESB(*args):
    tool_pred = args[0]
    dna_regions = args[1]
    FGENESB_ORFs = collections.OrderedDict()
    for dna_region in dna_regions:
        FGENESB_ORFs[dna_region] = collections.OrderedDict()
    for dna_region in dna_regions:
        genome = dna_regions[dna_region][0]
        genome_size = len(genome)
        genome_rev = revCompIterative(genome)
        with open(tool_pred, 'r') as FGENESB_input:
            for line in FGENESB_input:
                if '>GENE' in line:
                    line = line.split()
                    if '2208' in line:
                        print("ss")
                    if len(line) == 10 and dna_region in line[0] and ">GENE" in line[0]:
                        start = int(line[2])
                        stop = int(line[4])
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
                        orf = [strand, startCodon, stopCodon, 'CDS', 'FGENESB']
                        FGENESB_ORFs.update({po: orf})

    for group in FGENESB_ORFs:
        FGENESB_ORFs[group] = sortORFs(FGENESB_ORFs[group])
    return FGENESB_ORFs
