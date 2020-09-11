import collections

from ..utils import revCompIterative


def MetaGeneAnnotator(input_to_analyse,Genome):
    MetaGeneAnnotator_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    with open('Tools/MetaGeneAnnotator/'+input_to_analyse,'r') as MetaGeneAnnotator_input:
        for line in MetaGeneAnnotator_input:
            line = line.split()
            if len(line) == 11:
                if "gene_" in line[0]:
                    start = int(line[1])
                    stop = int(line[2])
                    strand = line[3]
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
                    MetaGeneAnnotator_ORFs.update({po:orf})
    return MetaGeneAnnotator_ORFs

