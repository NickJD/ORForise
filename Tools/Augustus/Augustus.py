import collections
from ..utils import revCompIterative

def Augustus(input_to_analyse,Genome):
    Augustus_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/Augustus/' + input_to_analyse, 'r') as Augustus_input:
        for line in Augustus_input:
            line = line.split()
            if len(line) == 12 and "CDS" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
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
                Augustus_ORFs.update({po: orf})
    return Augustus_ORFs

































