import collections

import sys

sys.path.append('../')
from ..utils import revCompIterative

def StORF(input_to_analyse,genome):
    storf_orfs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)
    with open('Tools/StORF/'+input_to_analyse,'r') as storf_input:
        for line in storf_input:
            if "#" not in line:
                line = line.split()
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
                if '-' in strand:
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
    return storf_orfs








