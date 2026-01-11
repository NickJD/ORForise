import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def GeneMark_HA(*args):
    tool_pred = args[0]
    dna_regions = args[1]
    if not dna_regions: # This triggers if dna_regions is an empty dict (GFF_Intersect passed nothing)
        dna_regions = collections.OrderedDict()
        with open(tool_pred, 'r') as GeneMarkHA_input:
            for line in GeneMarkHA_input:
                line = line.split()
                if len(line) >= 9 and "CDS" in line[5] and line[0] not in dna_regions:
                    dna_regions[line[0]] = []  # Placeholder for genome sequence
        return dna_regions

    geneMarkHA_ORFs = collections.OrderedDict()
    for dna_region in dna_regions:
        geneMarkHA_ORFs[dna_region] = collections.OrderedDict()
    for dna_region in dna_regions:
        try:
            genome = dna_regions[dna_region][0]
        except IndexError:
            genome = dna_regions[dna_region]
        genome_size = len(genome)
        genome_rev = revCompIterative(genome)
        with open(tool_pred, 'r') as GeneMarkHA_input:
            for line in GeneMarkHA_input:
                line = line.split()
                if len(line) >= 9 and "CDS" in line[5] and dna_region in line[0]:
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
                    orf = [strand, startCodon, stopCodon, 'CDS', 'GeneMarkHA']
                    geneMarkHA_ORFs.update({po: orf})

    for group in geneMarkHA_ORFs:
        geneMarkHA_ORFs[group] = sortORFs(geneMarkHA_ORFs[group])
    return geneMarkHA_ORFs
