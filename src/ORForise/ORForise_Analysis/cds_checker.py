from importlib import import_module
import argparse
import collections
from datetime import date
import sys
try:
    from ORForise.src.ORForise.utils import revCompIterative
except ImportError:

     from ORForise.utils import revCompIterative



parser = argparse.ArgumentParser()
parser.add_argument('-dna', '--genome_dna', required=True, help='Genome DNA file (.fa) which both annotations '
                                                                'are based on')
parser.add_argument('-gff', '--genome_gff', required=True,
                    help='Which annotation file to add to reference annotation?')
args = parser.parse_args()




def cds_checker(genome_dna,genome_gff):
    genome_seq = ""
    with open(genome_dna, 'r') as genome_fasta:
        for line in genome_fasta:
            line = line.replace("\n", "")
            if not line.startswith('>'):
                genome_seq += str(line)
            else:
                genome_id = line.split()[0].replace('>','')

    ###########################################
    genome_size = len(genome_seq)
    genome_rev = revCompIterative(genome_seq)
    cds_dict = collections.OrderedDict()  # Order is important
    count = 0
    with open(genome_gff, 'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "biological_region" in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]

                    if '-' in strand:  # Reverse Compliment starts and stops adjusted
                        r_start = genome_size - stop
                        r_stop = genome_size - start
                        startCodon = genome_rev[r_start:r_start + 3]
                        stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                        length = abs(start - stop-1)
                    elif '+' in strand:
                        startCodon = genome_seq[start - 1:start + 2]
                        stopCodon = genome_seq[stop - 3:stop]
                        length = abs(start-1 - stop)
                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon]
                    cds_dict.update({po: orf})

                    if length % 3 == 0:
                        print("In-Fame")
                    else:
                        sys.exit("W")


                elif "bio" in line[2]:
                    sys.exit("SAS")
            except IndexError:
                continue


if __name__ == "__main__":
    cds_checker(**vars(args))

    print("Complete")
