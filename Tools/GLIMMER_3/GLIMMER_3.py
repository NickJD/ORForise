import collections
import sys
sys.path.append('../')
from DNA_Reverse_Compliment import revCompIterative




def GLIMMER_3(input_to_analyse,Genome):
    GLIMMER_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    #GLIMMER reverses the start and stop positions for ORFS on the 3-prime-end
    with open('GLIMMER_3/' + input_to_analyse, 'r') as glimmer_input:
        for line in glimmer_input:
            if '>' not in line:
                line = line.split()
                if "orf" in line[0]:
                    if '-' in line[3]:  # Reverse Compliment starts and stops to confirm to our definition
                        # Switched to match Sense Strand
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
                        startCodon = Genome[start - 1:start - 1 + 3]
                        stopCodon = Genome[stop - 3:stop - 1 + 1]

                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon]
                    GLIMMER_ORFs.update({po: orf})

        # TO make sure all ORFs are in order - Not WOrking
    #GLIMMER_ORFs = collections.OrderedDict(sorted(GLIMMER_ORFs.items(), key=lambda t: t[0]))
    return GLIMMER_ORFs

