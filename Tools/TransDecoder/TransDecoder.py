import collections
from ..DNA_Reverse_Compliment import revCompIterative

def TransDecoder(input_to_analyse,Genome):
    transDecoder_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/TransDecoder/'+input_to_analyse,'r') as transDecoder_Input:
        for line in transDecoder_Input:
            line = line.split()
            if len(line) == 9 and "transdecoder" in line[1] and "CDS" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                if '-' in strand: # Reverse Compliment starts and stops adjusted
                    r_start = Genome_Size - stop
                    r_stop = Genome_Size - start
                    startCodon = Genome_rev[r_start:r_start + 3]
                    stopCodon = Genome_rev[r_stop - 2:r_stop + 1]
                elif '+' in strand:
                    startCodon = Genome[start - 1:start +2]
                    stopCodon = Genome[stop - 3:stop]
                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                transDecoder_ORFs.update({po:orf})
    return transDecoder_ORFs



















