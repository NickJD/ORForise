import collections
import sys
sys.path.append('../Comparison')
from DNA_Reverse_Compliment import revCompIterative

def Augustus(input_to_analyse,Genome):

    Augustus_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)

    with open('Augustus/' + input_to_analyse, 'r') as Augustus_input:
        for line in Augustus_input:
            if "Chromosome	AUGUSTUS	CDS" in line:
                line = line.split()
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
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
                Augustus_ORFs.update({po: orf})

    # TO make sure all ORFs are in order - Not Working
    #Augustus_ORFs = collections.OrderedDict(sorted(Augustus_ORFs.items(), key=lambda t: t[0]))
    return Augustus_ORFs

































