import collections

from ..utils import revCompIterative
from ..utils import sortORFs

def FGENESB(input_to_analyse,Genome):
    FGENESB_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/FGENESB/' + input_to_analyse, 'r') as FGENESB_input:
        for line in FGENESB_input:
            if '>GENE' in line:
                line = line.split()
                if '2208' in line:
                    print("ss")
                if len(line) == 10 and ">GENE" in line[0]:
                    start = int(line[2])
                    stop = int(line[4])
                    strand = line[9]
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

    FGENESB_ORFs = sortORFs(FGENESB_ORFs)
    return FGENESB_ORFs

































