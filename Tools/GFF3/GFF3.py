import collections

from ..utils import revCompIterative
from ..utils import sortORFs

def GFF3(genome_to_compare,parameters,genome):
    GFF3_ORFs = collections.OrderedDict()
    genome_Size = len(genome)
    genome_rev = revCompIterative(genome)
    with open('Tools/GFF3/GFF3_'+genome_to_compare+'.gff','r') as prodigal_input:
        for line in prodigal_input:
            if '#' not in line:
                line = line.split()
                if "CDS" in line[2]: #and len(line) ==9?
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
                    GFF3_ORFs.update({po:orf})

    GFF3_ORFs = sortORFs(GFF3_ORFs)
    return GFF3_ORFs