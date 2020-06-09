import collections
import sys
sys.path.append('../')
from DNA_Reverse_Compliment import revCompIterative




def TransDecoder(input_to_analyse,Genome):

    transDecoder_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)


    with open('TransDecoder/'+input_to_analyse,'r') as transDecoder_Input:
        for line in transDecoder_Input:
            if "Chromosome	transdecoder	CDS" in line :
                start = int(line.split('\t')[3])
                stop = int(line.split('\t')[4])
                strand = line.split()[6]
                if '-' in line.split('\t')[6]: #Reverse Compliment starts and stops to confirm to our deffinition
                    #Switched to match Sense Strand
                    r_start = Genome_Size - stop
                    r_stop = Genome_Size - start
                    startCodon = Genome_rev[r_start:r_start + 3]
                    stopCodon = Genome_rev[r_stop - 2:r_stop + 1]
                elif '+' in line.split('\t')[6]:
                    startCodon = Genome[start - 1:start -1 + 3]
                    stopCodon = Genome[stop - 3:stop -1 + 1]


                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                transDecoder_ORFs.update({po:orf})

    #TO make sure all ORFs are in order - Not working
    #ransDecoder_ORFs = collections.OrderedDict(sorted(transDecoder_ORFs.items(),key=lambda t: t[0]))
    return transDecoder_ORFs



















