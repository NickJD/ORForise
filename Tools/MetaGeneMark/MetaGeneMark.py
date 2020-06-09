import collections
import sys
sys.path.append('../Comparison')
from DNA_Reverse_Compliment import revCompIterative


def MetaGeneMark(input_to_analyse,Genome):

    metaGeneMarkORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)

    with open('MetaGeneMark/'+input_to_analyse,'r') as metaGeneMark_input:
        for line in metaGeneMark_input:
            if "Chromosome" and "CDS" in line:
                start = int(line.split('\t')[3])
                stop = int(line.split('\t')[4])
                strand = line.split()[9]
                if '-' in line.split('\t')[6]: #Reverse Compliment starts and stops to confirm to our deffinition
                    #Switched to match Sense Strand
                    r_start = Genome_Size - stop
                    r_stop = Genome_Size - start
                    startCodon = Genome_rev[r_start:r_start + 3]
                    stopCodon = Genome_rev[r_stop - 2:r_stop + 1]
                elif '+' in line.split('\t')[6]:
                    startCodon = Genome[start - 1:start - 1 + 3]
                    stopCodon = Genome[stop - 3:stop - 1 + 1]

                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                metaGeneMarkORFs.update({po:orf})







    #Sort not working
    #metaGeneMarkORFs  = collections.OrderedDict(sorted(metaGeneMarkORFs .items(),key=lambda t: t[0]))
    return metaGeneMarkORFs




