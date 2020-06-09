import collections
import sys
sys.path.append('../')
from DNA_Reverse_Compliment import revCompIterative


def FragGeneScan(input_to_analyse, Genome):


    FragGeneScan_ORFs = collections.OrderedDict()
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)

    with open('FragGeneScan/'+input_to_analyse,'r') as fragGeneScan_input:
        for line in fragGeneScan_input:
            if "Chromosome	FGS	CDS" in line :
                start = int(line.split()[3])
                stop = int(line.split()[4])
                strand = line.split()[6]
                if '-' in line.split()[6]:  # Reverse Compliment starts and stops to confirm to our definition
                    # Switched to match Sense Strand
                    r_start = Genome_Size - stop
                    r_stop = Genome_Size - start
                    startCodon = Genome_rev[r_start:r_start + 3]
                    stopCodon = Genome_rev[r_stop - 2:r_stop + 1]


                elif '+' in line.split()[6]:
                    startCodon = Genome[start - 1:start - 1 + 3]
                    stopCodon = Genome[stop - 3:stop - 1 + 1]

                po = str(start) + ',' + str(stop)
                orf = [strand, startCodon, stopCodon]
                FragGeneScan_ORFs.update({po: orf})

                # TO make sure all ORFs are in order - Not Working
    #FragGeneScan_ORFs = collections.OrderedDict(sorted(FragGeneScan_ORFs.items(), key=lambda t: t[0]))
    return FragGeneScan_ORFs