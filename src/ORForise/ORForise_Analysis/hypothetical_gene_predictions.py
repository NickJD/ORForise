import argparse
import collections



parser = argparse.ArgumentParser()
parser.add_argument('-ref', '--reference_annotation', required=True,
                    help='Which reference annotation file to use as reference?')
parser.add_argument('-tp', '--tool_prediction', required=True, help='Tool genome prediction file')
args = parser.parse_args()


def main(reference_annotation,tool_prediction):
    ref_genes = collections.OrderedDict()  # Order is important
    num_hypo = 0
    with open(reference_annotation, 'r') as genome_gff:
        for line in genome_gff:
            line = line.split('\t')
            try:
                if "gene" in line[2] and len(line) == 9: # Have to use gene and not CDS here because the CDS tag does not contain the classification
                    start = line[3]
                    stop = line[4]
                    gene = start+'_'+stop
                    if "hypothetical" in line[8]:
                        ref_genes.update({gene:"hypothetical"})
                        num_hypo +=1
                    else:
                        ref_genes.update({gene: "N/A"})
            except IndexError:
                continue
    print("Number of hypothetic genes: "+str(num_hypo))
    ###############################
    perfect_match_hypo, partial_match_hypo, missed_hypo = 0,0,0
    perfect_match = False
    partial_match = False
    missed = False
    posss = []
    with open(tool_prediction, 'r') as tool_in:
        for line in tool_in:
            if line.startswith("Perfect_Match_Genes:"):
                perfect_match = True
            elif line.startswith("Partial_Match_Genes:"):
                perfect_match = False
                partial_match = True
            elif line.startswith("Missed_Genes:"):
                perfect_match = False
                partial_match = False
                missed = True
            elif line.startswith("ORF_Without_Corresponding_Gene_in_Ensembl"):
                break
            #################
            if perfect_match == True:
                if line.startswith('>'):
                    pos = line.split('_')[1] +'_'+ line.split('_')[2]
                    try:
                        if "hypothetical" in ref_genes[pos]:
                            perfect_match_hypo +=1
                            print(pos)
                            if pos in posss:
                                print("WE")
                            posss.append(pos)
                    except KeyError:
                        continue
            elif partial_match == True:
                if line.startswith('Gene:'): # Different tags
                    pos = line.split('_')[0].split(':')[1] +'_'+ line.split('_')[1] # should change Orforise output
                    try:
                        if "hypothetical" in ref_genes[pos]:
                            partial_match_hypo +=1
                            print(pos)
                    except KeyError:
                        continue
            elif missed == True:
                if line.startswith('>'):
                    pos = line.split('_')[1] +'_'+ line.split('_')[2]
                    try:
                        if "hypothetical" in ref_genes[pos]:
                            missed_hypo +=1
                            print(pos)
                    except KeyError:
                        continue

    print("finished")

if __name__ == "__main__":
    main(**vars(args))

    print("Complete")
