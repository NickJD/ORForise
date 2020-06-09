import collections
import sys
sys.path.append('../Comparison')
from DNA_Reverse_Compliment import revCompIterative

def GeneMark_HMM(input_to_analyse,Genome):

    GeneMark_HMM_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)

    with open('GeneMark_HMM/'+input_to_analyse,'r') as GeneMark_HMM_input:
            for line in GeneMark_HMM_input:
                line = line.split()
                if len(line) >= 9 and "Chromosome" in line[0] and "CDS" in line[5]:
                    start = int(line[6])
                    stop = int(line[7])
                    strand = line[9]
                    if '-' in strand:  # Reverse Compliment starts and stops to confirm to our definition
                        # Switched to match Sense Strand
                        r_start = Genome_Size - stop
                        r_stop = Genome_Size - start
                        startCodon = Genome_rev[r_start:r_start + 3]
                        stopCodon = Genome_rev[r_stop - 2:r_stop + 1]


                    elif '+' in strand:
                        startCodon = Genome[start - 1:start - 1 + 3]
                        stopCodon = Genome[stop - 3:stop - 1 + 1]

                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon]
                    GeneMark_HMM_ORFs.update({po: orf})

                    # TO make sure all ORFs are in order - Not Working
    #GeneMark_HMM_ORFs = collections.OrderedDict(sorted(GeneMark_HMM_ORFs.items(), key=lambda t: t[0]))
    return GeneMark_HMM_ORFs

















