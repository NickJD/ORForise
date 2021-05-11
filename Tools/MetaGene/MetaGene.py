import collections

from ORForise.utils import revCompIterative
from ORForise.utils import sortORFs

def MetaGene(genome_to_compare,parameters,genome):
    metaGene_ORFs = collections.OrderedDict()
    genome_Size = len(genome)
    genome_rev = revCompIterative(genome)
    with open('Tools/MetaGene/MetaGene_'+genome_to_compare+'.txt','r') as MetaGene_input:
        for line in MetaGene_input:
            line = line.split()
            if len(line) >= 6 and ("-" in line or '+' in line):
                start = int(line[0])
                stop = int(line[1])
                strand = line[2]
                if '-' in strand:  # Reverse Compliment starts and stops adjusted
                    r_start = genome_Size - stop
                    r_stop = genome_Size - start
                    startCodon = genome_rev[r_start:r_start + 3]
                    stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                elif '+' in strand:
                    startCodon = genome[start - 1:start+2]
                    stopCodon = genome[stop - 3:stop]
                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                metaGene_ORFs.update({po:orf})

    metaGene_ORFs = sortORFs(metaGene_ORFs)
    return metaGene_ORFs

