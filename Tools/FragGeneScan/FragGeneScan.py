import collections

from ..utils import revCompIterative
from ..utils import sortORFs

def FragGeneScan(genome_to_compare,parameters,genome):
    fragGeneScan_ORFs = collections.OrderedDict()
    genome_Size = len(genome)
    genome_rev = revCompIterative(genome)
    with open('Tools/FragGeneScan/FragGeneScan_'+genome_to_compare+'.gff','r') as fragGeneScan_input:
        for line in fragGeneScan_input:
            line = line.split()
            if len(line) == 10 and "FGS" in line[1] and "CDS" in line [2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                if '-' in strand:  # Reverse Compliment starts and stops adjusted
                    r_start = genome_Size - stop
                    r_stop = genome_Size - start
                    startCodon = genome_rev[r_start:r_start + 3]
                    stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                elif '+' in strand:
                    startCodon = genome[start - 1:start+2]
                    stopCodon = genome[stop - 3:stop]
                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                fragGeneScan_ORFs.update({po: orf})

    fragGeneScan_ORFs = sortORFs(fragGeneScan_ORFs)
    return fragGeneScan_ORFs