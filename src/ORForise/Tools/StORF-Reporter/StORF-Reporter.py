import collections

try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def StORF_Reporter(*args):
    tool_pred = args[0]
    dna_regions = args[1]
    if not dna_regions: # This triggers if dna_regions is an empty dict (GFF_Intersect passed nothing)
        dna_regions = collections.OrderedDict()
        with open(tool_pred, 'r') as StORF_Reporter_input:
            for line in StORF_Reporter_input:
                line = line.split()
                if 'StORF-Reporter' in line[1] or 'StoRF_Reporter' in line[1]  or 'StORF' in line[1] or 'StORF-Reporter' in line[1] and line[0] not in dna_regions:
                    dna_regions[line[0]] = []  # Placeholder for genome sequence
        return dna_regions

    storf_ORFs = collections.OrderedDict()
    for dna_region in dna_regions:
        storf_ORFs[dna_region] = collections.OrderedDict()
    for dna_region in dna_regions:
        try:
            genome = dna_regions[dna_region][0]
        except IndexError:
            genome = dna_regions[dna_region]
        genome_size = len(genome)
        genome_rev = revCompIterative(genome)
        with open(tool_pred, 'r') as storf_input:
            for line in storf_input:
                if not line.startswith('#') and not line.startswith('\n'):
                    line = line.split()
                    if 'StORF-Reporter' in line[1] or 'StoRF_Reporter' in line[1]  or 'StORF' in line[1] or 'StORF-Reporter' in line[1] and dna_region in line[0]: # need to harmonise this.
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
                            startCodon = genome[start:start + 3]
                            stopCodon = genome[stop - 3:stop]
                        po = str(start) + ',' + str(stop)
                        orf = [strand, startCodon, stopCodon, 'CDS', 'StORF-Reporter'] # StORF/Con-StORF or CDS??
                        storf_ORFs.update({po: orf})

    for group in storf_ORFs:
        storf_ORFs[group] = sortORFs(storf_ORFs[group])
    return storf_ORFs
