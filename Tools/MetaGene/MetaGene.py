import collections
import sys
sys.path.append('../')
from DNA_Reverse_Compliment import revCompIterative

def MetaGene(input_to_analyse,Genome):
    MetaGene_ORFs = collections.OrderedDict()

    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)

    with open('MetaGene/'+input_to_analyse,'r') as MetaGene_input:
        for line in MetaGene_input:
            line = line.split()
            if "-" in line or '+' in line:
                start = int(line[0])
                stop = int(line[1])
                strand = line[2]
                if '-' in strand:  # Reverse Compliment starts and stops to confirm to our definition
                    # Switched to match Sense Strand
                    r_start = Genome_Size - stop
                    r_stop = Genome_Size - start
                    startCodon = Genome_rev[r_start:r_start + 3]
                    stopCodon = Genome_rev[r_stop - 2:r_stop + 1]


                elif '+' in strand:
                    startCodon = Genome[start - 1:start -1 + 3]
                    stopCodon = Genome[stop - 3:stop -1 + 1]

                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                MetaGene_ORFs.update({po:orf})

    #TO make sure all ORFs are in order - Not Working
    #MetaGene_ORFs = collections.OrderedDict(sorted(MetaGene_ORFs.items(),key=lambda t: t[0]))
    return MetaGene_ORFs

