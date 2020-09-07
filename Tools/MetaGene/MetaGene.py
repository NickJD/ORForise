import collections

from ..utils import revCompIterative

def MetaGene(input_to_analyse,Genome):
    MetaGene_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/MetaGene/'+input_to_analyse,'r') as MetaGene_input:
        for line in MetaGene_input:
            line = line.split()
            if len(line) >= 6 and ("-" in line or '+' in line):
                start = int(line[0])
                stop = int(line[1])
                strand = line[2]
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
                MetaGene_ORFs.update({po:orf})
    return MetaGene_ORFs

