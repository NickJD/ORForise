import collections

from ..utils import revCompIterative
from ..utils import sortORFs

def GeneMark_HMM(genome_to_compare,parameters,genome):
    geneMark_HMM_ORFs = collections.OrderedDict()
    genome_Size = len(genome)
    genome_rev = revCompIterative(genome)
    with open('Tools/GeneMark_HMM/GeneMark_HMM_'+genome_to_compare+'_'+parameters+'.gff','r') as GeneMark_HMM_input:
            for line in GeneMark_HMM_input:
                line = line.split()
                if len(line) >= 9 and "CDS" in line[5]:
                    start = int(line[6])
                    stop = int(line[7])
                    strand = line[9]
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
                    geneMark_HMM_ORFs.update({po: orf})

    geneMark_HMM_ORFs = sortORFs(geneMark_HMM_ORFs)
    return geneMark_HMM_ORFs

















