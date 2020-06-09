import collections
import sys
sys.path.append('../Comparison')
from DNA_Reverse_Compliment import revCompIterative


def GeneMark_HA(input_to_analyse,Genome):

    GeneMark_HA_ORFs = collections.OrderedDict()

    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)

    with open('GeneMark_HA/'+input_to_analyse,'r') as GeneMark_HA_input:
        for line in GeneMark_HA_input:
            line = line.split()
            if "Chromosome" in line:
                start = int(line[6])
                stop = int(line[7])
                strand = line[9]
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
                GeneMark_HA_ORFs.update({po: orf})

    #TO make sure all ORFs are in order - Not Working
    #GeneMark_HA_ORFs = collections.OrderedDict(sorted(GeneMark_HA_ORFs.items(), key=lambda t: t[0]))
    return GeneMark_HA_ORFs





