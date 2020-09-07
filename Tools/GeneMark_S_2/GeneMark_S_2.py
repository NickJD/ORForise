import collections
from ..utils import revCompIterative

def GeneMark_S_2(input_to_analyse,Genome):
    GeneMark_S_2_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/GeneMark_S_2/'+input_to_analyse,'r') as GeneMark_S_2_input:
        for line in GeneMark_S_2_input:
            line = line.split()
            if len(line) >= 9 and "CDS" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                if '-' in strand:  # Reverse Compliment starts and stops adjusted
                    r_start = Genome_Size - stop
                    r_stop = Genome_Size - start
                    startCodon = Genome_rev[r_start:r_start + 3]
                    stopCodon = Genome_rev[r_stop - 2:r_stop + 1]
                elif '+' in strand:
                    startCodon = Genome[start - 1:start + 2]
                    stopCodon = Genome[stop - 3:stop]
                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                GeneMark_S_2_ORFs.update({po:orf})
    return GeneMark_S_2_ORFs


