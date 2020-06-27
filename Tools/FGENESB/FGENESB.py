import collections
from ..DNA_Reverse_Compliment import revCompIterative

def FGENESB(input_to_analyse,Genome):
    FGENESB_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/FGENESB/' + input_to_analyse, 'r') as FGENESB_input:
        for line in FGENESB_input:
            line = line.split()
            if len(line) == 11 and "CDS" in line[6]:
                start = int(line[7])
                stop = int(line[9])
                strand = line[5]
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
                FGENESB_ORFs.update({po: orf})
            elif 'Predicted protein(s):' in line:
                break
    return FGENESB_ORFs

































