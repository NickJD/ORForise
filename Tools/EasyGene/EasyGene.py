import collections
from ..DNA_Reverse_Compliment import revCompIterative

def EasyGene(input_to_analyse,Genome):
    EasyGene_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/EasyGene/' + input_to_analyse, 'r') as EasyGene_input:
        for line in EasyGene_input:
            line = line.split()
            if len(line) == 10  and "CDS" in line[2]:
                start = int(line[3])
                stop = int(line[4])
                strand = line[6]
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
                EasyGene_ORFs.update({po: orf})
    return EasyGene_ORFs

