import collections
from ..DNA_Reverse_Compliment import revCompIterative

def GLIMMER_3(input_to_analyse,Genome):
    GLIMMER_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/GLIMMER_3/' + input_to_analyse, 'r') as glimmer_input: #GLIMMER_3 reverses the start and stop positions for ORFS on the negative strand
        for line in glimmer_input:
            if '>' not in line: # This will not work with multiple contigs
                line = line.split()
                if len(line) == 5 and "orf" in line[0]:
                    if '-' in line[3]:  # Reverse Compliment starts and stops adjusted -  Switched to match Sense Strand
                        start = int(line[2])
                        stop = int(line[1])
                        strand = '-'
                        r_start = Genome_Size - stop
                        r_stop = Genome_Size - start
                        startCodon = Genome_rev[r_start:r_start + 3]
                        stopCodon = Genome_rev[r_stop - 2:r_stop + 1]
                    elif '+' in line[3]:
                        start = int(line[1])
                        stop = int(line[2])
                        strand = '+'
                        startCodon = Genome[start - 1:start+3]
                        stopCodon = Genome[stop - 3:stop]
                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon]
                    GLIMMER_ORFs.update({po: orf})
    return GLIMMER_ORFs

