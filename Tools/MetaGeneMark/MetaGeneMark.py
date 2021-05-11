import collections

from ORForise.utils import revCompIterative
from ORForise.utils import sortORFs

def MetaGeneMark(genome_to_compare,parameters,genome):
    metaGeneMarkORFs = collections.OrderedDict()
    genome_Size = len(genome)
    genome_rev = revCompIterative(genome)
    with open('Tools/MetaGeneMark/MetaGeneMark_'+genome_to_compare+'.gff','r') as metaGeneMark_input:
        for line in metaGeneMark_input:
            line = line.split()
            if len(line) == 19:
                if 'GeneMark.hmm' in line[4] and "CDS" in line[5]:
                    start = int(line[6])
                    stop = int(line[7])
                    strand = line[9]
                    if '-' in strand: # Reverse Compliment starts and stops adjusted
                        r_start = genome_Size - stop
                        r_stop = genome_Size - start
                        startCodon = genome_rev[r_start:r_start + 3]
                        stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                    elif '+' in strand:
                        startCodon = genome[start - 1:start+2]
                        stopCodon = genome[stop - 3:stop]
                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon]
                    metaGeneMarkORFs.update({po:orf})

    metaGeneMarkORFs = sortORFs(metaGeneMarkORFs)
    return metaGeneMarkORFs
