import collections
from ..utils import revCompIterative

def GeneMark_HMM(input_to_analyse,Genome):
    GeneMark_HMM_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/GeneMark_HMM/'+input_to_analyse,'r') as GeneMark_HMM_input:
            for line in GeneMark_HMM_input:
                line = line.split()
                if len(line) >= 9 and "CDS" in line[5]:
                    start = int(line[6])
                    stop = int(line[7])
                    strand = line[9]
                    if '-' in strand:  # Reverse Compliment starts and stops adjusted
                        r_start = Genome_Size - stop
                        r_stop = Genome_Size - start
                        startCodon = Genome_rev[r_start:r_start + 3]
                        stopCodon = Genome_rev[r_stop - 2:r_stop + 1]
                    elif '+' in strand:
                        startCodon = Genome[start - 1:start+2]
                        stopCodon = Genome[stop - 3:stop]
                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon]
                    GeneMark_HMM_ORFs.update({po: orf})
    return GeneMark_HMM_ORFs

















