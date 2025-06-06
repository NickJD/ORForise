import collections
import sys
try:
    from utils import revCompIterative
    from utils import sortORFs
except ImportError:
    from ORForise.utils import revCompIterative
    from ORForise.utils import sortORFs


def GFF(*args):
    tool_pred = args[0]
    genome = args[1]
    #types = args[2]
    GFF_ORFs = collections.OrderedDict()
    genome_size = len(genome)
    genome_rev = revCompIterative(genome)
    with open(tool_pred, 'r') as gff_input:
        for line in gff_input:
            if '#' not in line:
                line = line.split('\t')
                #gene_types = types.split(',') - Temporary fix
                #if any(gene_type == line[2] for gene_type in gene_types) and len(line) == 9:  # line[2] for normalrun
                if 'CDS' in line[2] and len(line) == 9:
                    start = int(line[3])
                    stop = int(line[4])
                    strand = line[6]
                    info = line[8]
                    if stop >= genome_size:
                        extra_stop = stop - genome_size
                        corrected_stop = genome_size
                        if '-' in strand:  # Reverse Compliment starts and stops adjusted
                            r_start = genome_size - corrected_stop
                            r_stop = genome_size - start
                            seq = genome_rev[r_start:r_stop + 1]
                            extra_seq = genome_rev[-extra_stop - 1:]
                            seq = extra_seq+seq
                            startCodon = seq[:3]
                            stopCodon = seq[-3:]
                        elif '+' in strand:
                            seq = genome[start -1 :corrected_stop]
                            extra_seq = genome[:extra_stop +1]
                            seq = seq+extra_seq
                            startCodon = seq[:3]
                            stopCodon = seq[-3:]
                    else:
                        if '-' in strand:  # Reverse Compliment starts and stops adjusted
                            r_start = genome_size - stop
                            r_stop = genome_size - start
                            startCodon = genome_rev[r_start:r_start + 3]
                            stopCodon = genome_rev[r_stop - 2:r_stop + 1]
                        elif '+' in strand:
                            startCodon = genome[start - 1:start + 2]
                            stopCodon = genome[stop - 3:stop]
                    po = str(start) + ',' + str(stop)
                    orf = [strand, startCodon, stopCodon, line[2],info] # This needs to detect the type
                    GFF_ORFs.update({po: orf})
                # elif "CDS" in line[2]:
                #     sys.exit("SAS")

    GFF_ORFs = sortORFs(GFF_ORFs)
    return GFF_ORFs
