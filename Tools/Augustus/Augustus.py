import collections

from ..utils import revCompIterative
from ..utils import sortORFs

def Augustus(genome_to_compare,parameters,genome):
    augustus_ORFs = collections.OrderedDict()
    genome_Size = len(genome)
    genome_rev = revCompIterative(genome)
    with open('Tools/Augustus/'+genome_to_compare+'_'+parameters+'.csv','r') as Augustus_input:
        for line in Augustus_input:
            line = line.split()
            if len(line) == 12 and "CDS" in line[2]:
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
                augustus_ORFs.update({po: orf})

    augustus_ORFs = sortORFs(augustus_ORFs)
    return augustus_ORFs

































