import collections
import sys
sys.path.append('../')
from ..DNA_Reverse_Compliment import revCompIterative



def STORF(input_to_analyse,genome):

    storf_orfs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)

    with open('Tools/STORF/'+input_to_analyse,'r') as storf_input:
        for line in storf_input:
            if "#" not in line:
                line = line.split()
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                if '-' in strand:  # Reverse Compliment starts and stops to confirm to our definition
                    # Switched to match Sense Strand
                    r_start = genome_size - stop
                    r_stop = genome_size - start
                    startcodon = genome_rev[r_start:r_start + 3]
                    stopcodon = genome_rev[r_stop - 2:r_stop + 1]

                elif '+' in strand:
                    startcodon = genome[start - 1:start -1 + 3]
                    stopcodon = genome[stop - 3:stop -1 + 1]

                po = str(start) + ',' + str(stop)
                orf = [strand, startcodon, stopcodon]
                storf_orfs.update({po:orf})

    #TO make sure all ORFs are in order
    #storf_orfs = collections.OrderedDict(sorted(storf_orfs.items(),key=lambda t: t[0]))
    return storf_orfs








