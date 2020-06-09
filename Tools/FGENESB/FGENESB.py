import collections
import sys
sys.path.append('../')
from DNA_Reverse_Compliment import revCompIterative

def FGENESB(input_to_analyse,Genome):

    FGENESB_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)

    with open('FGENESB/' + input_to_analyse, 'r') as FGENESB_input:
        for line in FGENESB_input:
            if "CDS" in line:
                line = line.split()
                start = int(line[7])
                stop = int(line[9])
                strand = line[5]
                if '-' in strand:  # Reverse Compliment starts and stops to confirm to our definition
                    # Switched to ma
                    # tch Sense Strand
                    r_start = Genome_Size - stop
                    r_stop = Genome_Size - start
                    startCodon = Genome_rev[r_start:r_start + 3]
                    stopCodon = Genome_rev[r_stop - 2:r_stop + 1]


                elif '+' in strand:
                    startCodon = Genome[start - 1:start - 1 + 3]
                    stopCodon = Genome[stop - 3:stop - 1 + 1]

                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                FGENESB_ORFs.update({po: orf})

            elif 'Predicted protein(s):' in line:
                break


            # TO make sure all ORFs are in order - Not Working
    #FGENESB_ORFs = collections.OrderedDict(sorted(FGENESB_ORFs.items(), key=lambda t: t[0]))
    return FGENESB_ORFs

































