import collections

from ORForise.utils import revCompIterative
from ORForise.utils import sortORFs

def GeneMark_S_2(genome_to_compare,parameters,genome):
    geneMark_S_2_ORFs = collections.OrderedDict()
    genome_Size = len(genome)
    genome_rev = revCompIterative(genome)
    with open('Tools/GeneMark_S_2/GeneMark_S_2_'+genome_to_compare+'.gff','r') as GeneMark_S_2_input:
        for line in GeneMark_S_2_input:
            line = line.split()
            if len(line) >= 9 and "CDS" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                if '-' in strand:  # Reverse Compliment starts and stops adjusted
                    r_start = genome_Size - stop
                    r_stop = genome_Size - start
                    startCodon = genome_rev[r_start:r_start + 3]
                    stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                elif '+' in strand:
                    startCodon = genome[start - 1:start + 2]
                    stopCodon = genome[stop - 3:stop]
                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                geneMark_S_2_ORFs.update({po:orf})

    geneMark_S_2_ORFs = sortORFs(geneMark_S_2_ORFs)
    return geneMark_S_2_ORFs


