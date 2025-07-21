import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def MetaGeneMark(*args):
    tool_pred = args[0]
    dna_regions = args[1]
    metaGeneMarkORFs = collections.OrderedDict()
    for dna_region in dna_regions:
        metaGeneMarkORFs[dna_region] = collections.OrderedDict()
    for dna_region in dna_regions:
        genome = dna_regions[dna_region][0]
        genome_size = len(genome)
        genome_rev = revCompIterative(genome)
        with open(tool_pred, 'r') as metaGeneMark_input:
            for line in metaGeneMark_input:
                line = line.split()
                if len(line) == 19:
                    if 'GeneMark.hmm' in line[4] and "CDS" in line[5] and dna_region in line[0]:
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
                        orf = [strand, startCodon, stopCodon, 'CDS', 'MetaGeneMark']
                        metaGeneMarkORFs.update({po: orf})

    for group in metaGeneMarkORFs:
        metaGeneMarkORFs[group] = sortORFs(metaGeneMarkORFs[group])
    return metaGeneMarkORFs
