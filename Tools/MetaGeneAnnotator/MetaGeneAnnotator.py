import collections
import sys
sys.path.append('../')
from DNA_Reverse_Compliment import revCompIterative

def MetaGeneAnnotator(input_to_analyse,Genome):
    MetaGeneAnnotator_ORFs = collections.OrderedDict()

    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)

    with open('MetaGeneAnnotator/'+input_to_analyse,'r') as MetaGeneAnnotator_input:
        for line in MetaGeneAnnotator_input:
            if "gene_" in line:
                line = line.split()
                start = int(line[1])
                stop = int(line[2])
                strand = line[3]
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
                MetaGeneAnnotator_ORFs.update({po:orf})

    #TO make sure all ORFs are in order - Not Working
    #MetaGeneAnnotator_ORFs = collections.OrderedDict(sorted(MetaGeneAnnotator_ORFs.items(),key=lambda t: t[0]))
    return MetaGeneAnnotator_ORFs

