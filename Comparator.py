import collections
import numpy as np
from Tools.DNA_Reverse_Compliment import  revCompIterative

class comparator:  # Class to hold global-type variables
    def __init__(self, perfect_Starts=0, perfect_Stops=0, perfect_Matches=0, genome_Seq='',
                 genome_Seq_Rev='',
                 genome_Size=0, correct_Frame_Number=0, expanded_Start=0,
                 expanded_Stop=0, expanded_CDS=0, matched_ORFs=collections.OrderedDict(),
                 unmatched_ORFs=collections.OrderedDict(), genes_Detected=collections.OrderedDict(),
                 genes_Undetected=collections.OrderedDict(),
                 out_Of_Frame_ORFs=collections.OrderedDict(), start_Difference=[], stop_Difference=[],
                 orf_Lengths=[], gene_Lengths=[], pos_Strand=0, neg_Strand=0):
        self.perfect_Starts, self.perfect_Stops, self.perfect_Matches, self.genome_Seq, self.genome_Seq_Rev, self.genome_Size, self.correct_Frame_Number, self.expanded_Start, self.expanded_Stop, self.expanded_CDS, \
        self.matched_ORFs, self.unmatched_ORFs, self.genes_Detected, self.genes_Undetected, self.out_Of_Frame_ORFs, self.start_Difference, self.stop_Difference, self.orf_Lengths, self.gene_Lengths, self.pos_Strand, \
        self.neg_Strand = perfect_Starts, perfect_Stops, perfect_Matches, genome_Seq, genome_Seq_Rev, genome_Size, correct_Frame_Number, expanded_Start, expanded_Stop, expanded_CDS, matched_ORFs, unmatched_ORFs, \
                          genes_Detected, genes_Undetected, out_Of_Frame_ORFs, start_Difference, stop_Difference, orf_Lengths, gene_Lengths, pos_Strand, neg_Strand

comp = comparator()


def orf_Unmatched(o_Start, o_Stop, o_Strand):
    if o_Strand == '-':
        r_Start = comp.genome_Size - o_Stop
        r_Stop = comp.genome_Size - o_Start
        Unmatched_ORF = str(o_Start) + ',' + str(o_Stop) + ',' + o_Strand + ',' +comp.genome_Seq_Rev[r_Start:r_Start+3] + ',' +comp.genome_Seq_Rev[r_Stop - 2:r_Stop + 1]
        seq = (comp.genome_Seq_Rev[r_Start-3:r_Stop+3])
        comp.unmatched_ORFs.update({Unmatched_ORF: seq})
    elif o_Strand == '+':
        Unmatched_ORF = str(o_Start) + ',' + str(o_Stop) + ',' + o_Strand + ',' + comp.genome_Seq[o_Start-1:o_Start+2] + ',' +comp.genome_Seq[o_Stop-3:o_Stop]
        seq = (comp.genome_Seq[o_Start - 1:o_Stop])
        comp.unmatched_ORFs.update({Unmatched_ORF: seq})

def genes_Unmatched(g_Start,g_Stop,g_Strand):
    if g_Strand == '-':
        r_Start = comp.genome_Size - g_Stop
        r_Stop = comp.genome_Size - g_Start
        missed_Gene = str(g_Start) + ',' + str(g_Stop) + ',' + g_Strand + ',' +comp.genome_Seq_Rev[r_Start:r_Start+3] + ',' +comp.genome_Seq_Rev[r_Stop - 2:r_Stop + 1]
        genSeq = (comp.genome_Seq_Rev[r_Start:r_Stop + 1])
        comp.genes_Undetected.update({missed_Gene: genSeq})
    elif g_Strand == '+':
        missed_Gene = str(g_Start) + ',' + str(g_Stop) + ',' + g_Strand + ',' + comp.genome_Seq[g_Start-1:g_Start+2] + ',' +comp.genome_Seq[g_Stop-3:g_Stop]
        genSeq = (comp.genome_Seq[g_Start - 1:g_Stop])
        comp.genes_Undetected.update({missed_Gene: genSeq})

def match_Statistics(o_Start,o_Stop,g_Start,g_Stop):
    if g_Start == o_Start:
        comp.perfect_Starts += 1
    if g_Stop == o_Stop:
        comp.perfect_Stops += 1
    ############ Calculate prediction precision
    comp.start_Difference.append(o_Start - g_Start)
    comp.stop_Difference.append(o_Stop - g_Stop)
    comp.correct_Frame_Number += 1
    if o_Start < g_Start and o_Stop > g_Stop:
        comp.expanded_CDS +=1
    if o_Start < g_Start:
        comp.expanded_Start +=1
    if o_Stop > g_Stop:
        comp.expanded_Stop +=1

def start_Codon_Count(orfs):
    atg,gtg,ttg,att,ctg,other = 0,0,0,0,0,0
    other_Starts = []
    for orf in orfs.values():
        codon = orf[1]
        if codon == 'ATG':
            atg += 1
        elif codon == 'GTG':
            gtg += 1
        elif codon == 'TTG':
            ttg += 1
        elif codon == 'ATT':
            att += 1
        elif codon == 'CTG':
            ctg += 1
        else:
            other += 1
            other_Starts.append(codon)
    atg_P = float((atg) * float(100) / float(len(orfs)))
    gtg_P = float((gtg) * float(100) / float(len(orfs)))
    ttg_P = float((ttg) * float(100) / float(len(orfs)))
    att_P = float((att) * float(100) / float(len(orfs)))
    ctg_P = float((ctg) * float(100) / float(len(orfs)))
    other_Start_P = float(other) * float(100) / float(len(orfs))
    return atg_P,gtg_P,ttg_P,att_P,ctg_P,other_Start_P,other_Starts

def stop_Codon_Count(orfs):
    tag,taa,tga,other = 0,0,0,0
    other_Stops = []
    for orf in orfs.values():
        codon = orf[2]
        if codon == 'TAG':
            tag += 1
        elif codon == 'TAA':
            taa += 1
        elif codon == 'TGA':
            tga += 1
        else:
            other += 1
            other_Stops.append(codon)
    tag_p = float((tag) * float(100) / float(len(orfs)))
    taa_p = float((taa) * float(100) / float(len(orfs)))
    tga_p = float((tga) * float(100) / float(len(orfs)))
    other_Stop_P = float(other) * float(100) / float(len(orfs))
    return tag_p,taa_p,tga_p,other_Stop_P,other_Stops

def candidate_ORF_Selection(gene_Set,candidate_ORFs): # Select ORF from candidates which is most similar to partially detected gene
    current_Coverage = 0
    candidate_ORF_Difference = 0
    pos = ''
    orf_Details = []
    for c_Pos, c_ORF_Details in candidate_ORFs.items():
        o_Start = int(c_Pos.split(',')[0])
        o_Stop = int(c_Pos.split(',')[1])
        coverage = c_ORF_Details[3]
        orf_Set = set(range(o_Start, o_Stop + 1))
        if coverage > current_Coverage:
            current_Coverage = coverage
            # Return set of elements outside the two sets/DNA ranges
            candidate_ORF_Difference = orf_Set.symmetric_difference(gene_Set)
            pos = c_Pos
            orf_Details = c_ORF_Details
        elif coverage == current_Coverage:
            current_ORF_Difference = orf_Set.symmetric_difference(gene_Set) # Pick least different ORF set from the Gene Set
            if len(current_ORF_Difference) > len(candidate_ORF_Difference):
                pos = c_Pos
                orf_Details = c_ORF_Details
        else:
            print("Match filtered out")
    return pos,orf_Details

def tool_comparison(genes,orfs,genome):
    comp.genome_Seq = genome
    comp.genome_Seq_Rev = revCompIterative(genome)
    comp.genome_Size = len(genome)
    #Loop through each gene to compare against predicted ORFs
    for gene_Num, gene_Details in genes.items():
        gene_Details = gene_Details.split(',')
        g_Start = int(gene_Details[0])
        g_Stop = int(gene_Details[1])
        g_Strand = gene_Details[2]
        gene_Length = g_Stop - g_Start
        comp.gene_Lengths.append(gene_Length)
        gene_Set = set(range(g_Start, g_Stop + 1))
        ####################### Checking Overlap of ORFs and pick best match - slow but confirms best match
        overlapping_ORFs = collections.OrderedDict()
        perfect_Match = False
        for pos, orf_Details  in orfs.items(): #Check if perfect match, if not check if match covers at least 75% of gene - Loop through ALL ORFs
            o_Start = int(pos.split(',')[0])
            o_Stop = int(pos.split(',')[1])
            o_Strand = orf_Details[0]
            orf_Set = set(range(o_Start, o_Stop + 1))
            if o_Stop <= g_Start or o_Start >= g_Stop: #Not caught up yet or ORFs may not be ordered
                continue
            elif o_Start == g_Start and o_Stop == g_Stop: #If perfect match, break and skip the rest of the ORFs
                perfect_Match = True
                break
            elif g_Start <= o_Start < g_Stop or g_Start < o_Stop < g_Stop: # If at least ORF Stop is within gene range
                overlap = len(gene_Set.intersection(orf_Set))
                coverage = 100 * float(overlap) / float(len(gene_Set))
                orf_Details.append(coverage)
                if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and coverage >= 75 :  #Only continue if ORF covers at least 75% of the gene and is in frame
                    overlapping_ORFs.update({pos:orf_Details})
                elif coverage >= 75:
                    comp.out_Of_Frame_ORFs.update({pos:orf_Details})
            elif o_Start <= g_Start and o_Stop >= g_Stop: # If ORF is extends one or both ends of the gene
                overlap = len(gene_Set.intersection(orf_Set))
                coverage = 100 * float(overlap) / float(len(gene_Set))
                orf_Details.append(coverage)
                if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and coverage >= 75:  #Only continue if ORF covers at least 75% of the gene and is in frame
                    overlapping_ORFs.update({pos:orf_Details})
                elif coverage >= 75:
                    comp.out_Of_Frame_ORFs.update({pos:orf_Details})
            else:
                print("Unexpected Error Finding ORFs")
        #Now Check that we select the best ORF
        #print(len(overlapping_ORFs))
        if perfect_Match == True: # Check if the ORF is a perfect match to the Gene
            comp.genes_Detected.update({str(gene_Details):pos})
            orf_Details.append(100)
            comp.matched_ORFs.update({pos:orf_Details})
            comp.perfect_Matches += 1
            comp.perfect_Starts += 1
            comp.perfect_Stops += 1
            comp.correct_Frame_Number += 1
            print('Perfect Match')
        elif perfect_Match == False and len(overlapping_ORFs) == 1: # If we do not have a perfect match but 1 ORF which has passed our filtering, we will calculate accordingly
            orf_Pos = list(overlapping_ORFs.keys())[0]
            o_Start = int(orf_Pos.split(',')[0])
            o_Stop = int(orf_Pos.split(',')[1])
            orf_Details = overlapping_ORFs[orf_Pos]
            comp.genes_Detected.update({str(gene_Details):orf_Pos})
            comp.matched_ORFs.update({orf_Pos: orf_Details})
            match_Statistics(o_Start,o_Stop,g_Start,g_Stop)
            print('Partial Match')
        elif perfect_Match == False and len(overlapping_ORFs) >= 1: # If we have more than 1 potential ORF match, we check to see which is the 'best' hit
            orf_Pos,orf_Details = candidate_ORF_Selection(gene_Set,overlapping_ORFs) # Return best match
            o_Start = int(orf_Pos.split(',')[0])
            o_Stop = int(orf_Pos.split(',')[1])
            comp.matched_ORFs.update({orf_Pos:orf_Details})
            comp.genes_Dected.update({str(gene_Details):orf_Pos})
            match_Statistics(o_Start,o_Stop,g_Start,g_Stop)
            print('There was more than 1 potential Match')
        elif len(comp.out_Of_Frame_ORFs) >=1: # Keep record of ORFs which overlap a gene but in the wrong frame
            print("Out of Frame ORF")
            genes_Unmatched(g_Start, g_Stop, g_Strand)
        else:
            genes_Unmatched(g_Start, g_Stop, g_Strand) # No hit
            print("No Hit")
    for key in comp.matched_ORFs: # Remove ORFs which are out of frame if ORF was correctly matched to another Gene
        if key in comp.out_Of_Frame_ORFs:
            del comp.out_Of_Frame_ORFs[key]
    print("Checked all predicted ORFs")
    print(comp.genes_Detected)
    min_Gene_Length = min(comp.gene_Lengths)
    max_Gene_Length = max(comp.gene_Lengths)
    median_Gene_Length = np.median(comp.gene_Lengths)
    ##################################################### ORF Lengths and Precision
    start_Difference = [x for x in comp.start_Difference if x != 0]
    stop_Difference = [x for x in comp.stop_Difference if x != 0]
    if len(start_Difference) >= 1:
        median_Start_Difference = np.median(start_Difference)
    else:
        median_Start_Difference = 0
    if len(stop_Difference) >= 1:
        median_Stop_Difference = np.median(stop_Difference)
    else:
        median_Stop_Difference = 0

    # Get Start and Stop Codon Usage
    atg_P, gtg_P, ttg_P, att_P, ctg_P, other_Start_P,other_Starts = start_Codon_Count(orfs)
    tag_P, taa_P, tga_P, other_Stop_P,other_Stops = stop_Codon_Count(orfs)
    # Count nucleotides found from ALL ORFs
    gene_Nuc_Count = np.zeros((comp.genome_Size), dtype=np.int)
    orf_Nuc_Count = np.zeros((comp.genome_Size), dtype=np.int)
    matched_ORF_Nuc_Count = np.zeros((comp.genome_Size), dtype=np.int)
    for g_Positions in genes.values():
        g_Start = int(g_Positions.split(',')[0])
        g_Stop = int(g_Positions.split(',')[1])
        gene_Nuc_Count[g_Start-1:g_Stop] = [1] # Changing all between the two positions to 1's
    for o_Positions,orf_Details in orfs.items():
        o_Start = int(o_Positions.split(',')[0])
        o_Stop = int(o_Positions.split(',')[1])
        o_Strand = orf_Details[0]
        comp.orf_Lengths.append(o_Stop - o_Start)
        orf_Nuc_Count[o_Start-1:o_Stop] = [1] # Changing all between the two positions to 1's
        # Get ORF Strand metrics:
        if o_Strand == "+":  # Get number of Positive and Negative strand ORFs
            comp.pos_Strand += 1
        elif o_Strand == "-":
            comp.neg_Strand += 1
        # Stats just for Unmatched ORFs
        if o_Positions not in list(comp.matched_ORFs.keys()):
            orf_Unmatched(o_Start, o_Stop, o_Strand)
    # Nucleotide Coverage calculated from ORFs matching a gene only
    for o_Positions,orf_Details in comp.matched_ORFs.items():
        o_Start = int(o_Positions.split(',')[0])
        o_Stop = int(o_Positions.split(',')[1])
        matched_ORF_Nuc_Count[o_Start-1 :o_Stop] = [1] # Changing all between the two positions to 1's
    gene_Coverage_Genome = 100 * float(np.count_nonzero(gene_Nuc_Count)) / float(comp.genome_Size)
    orf_Coverage_Genome = 100 * float(np.count_nonzero(orf_Nuc_Count)) / float(comp.genome_Size)
    matched_ORF_Coverage_Genome = 100 * float(np.count_nonzero(matched_ORF_Nuc_Count)) / float(comp.genome_Size)
    # gene and orf nucleotide union
    gene_ORF_Nuc_Union_Count = np.count_nonzero(gene_Nuc_Count & orf_Nuc_Count)
    #not gene but orf nucleotides
    not_Gene_Nuc = np.logical_not(gene_Nuc_Count) + [0 for i in range(len(gene_Nuc_Count))]
    not_Gene_Nuc_And_ORF_Count = np.count_nonzero(not_Gene_Nuc & orf_Nuc_Count)
    #not orf nucleotides but gene
    not_ORF_Nuc = np.logical_not(orf_Nuc_Count) + [0 for i in range(len(orf_Nuc_Count))]
    not_ORF_Nuc_And_Gene_Count = np.count_nonzero(not_ORF_Nuc & gene_Nuc_Count)
    #not gene or orf nucleotides
    not_Gene_Nuc_Not_ORF_Nuc_Count = np.count_nonzero(not_Gene_Nuc & not_ORF_Nuc)
    #Nucleotide 'accuracy'
    NT_TP = (float(gene_ORF_Nuc_Union_Count)  / float(np.count_nonzero(gene_Nuc_Count)))
    NT_FP = (float(not_Gene_Nuc_And_ORF_Count) / float(np.count_nonzero(gene_Nuc_Count)))
    NT_FN = (float(not_ORF_Nuc_And_Gene_Count) / float(np.count_nonzero(gene_Nuc_Count)))
    NT_TN = (float(not_Gene_Nuc_Not_ORF_Nuc_Count) / float(np.count_nonzero(gene_Nuc_Count)))
    NT_Precision = NT_TP / (NT_TP + NT_FP)
    NT_Recall = NT_TP / (NT_TP + NT_FN)
    NT_False_Discovery_Rate = NT_FP / (NT_FP + NT_TP)
    ################################# Precision and Recall of filtered ORFs
    TP = (len(comp.genes_Detected)  / len(genes))
    FP = (len(comp.unmatched_ORFs) / len(genes))
    FN = (len(comp.genes_Undetected)  / len(genes))
    try: # Incase no ORFs found a gene
        precision = TP/(TP+FP)
        recall = TP/(TP+FN)
        false_Discovery_Rate = FP/(FP+TP)
    except ZeroDivisionError:
        precision = 0
        recall = 0
        false_Discovery_Rate = 0
    min_ORF_Length = min(comp.orf_Lengths)
    max_ORF_Length = max(comp.orf_Lengths)
    median_ORF_Length = np.median(comp.orf_Lengths)

    # Percenting Metrics
    ORFs_Diff = (float(len(orfs)) - float(len(genes))) / float(len(genes)) * 100
    genes_Detected_Percentage = len(comp.genes_Detected) / len(genes) * 100
    matched_ORF_Percentage = len(comp.matched_ORFs) / len(orfs) * 100
    median_Length_Diff = ((float(median_ORF_Length) - median_Gene_Length) / median_Gene_Length) * 100
    min_Length_Diff = float((min_ORF_Length) - float(min_Gene_Length)) / float(min_Gene_Length) * 100
    max_Length_Diff = float((max_ORF_Length) - float(max_Gene_Length)) / float(max_Gene_Length) * 100
    pos_Strand_Percentage = (float(comp.pos_Strand) * 100 / float(len(orfs)))
    neg_Strand_Percentage = (float(comp.neg_Strand) * 100 / float(len(orfs)))
    #############################################
    try: # Incase no ORFs found a gene
        correct_Frame_Percentage = comp.correct_Frame_Number * 100 / (len(comp.genes_Detected) + len(comp.out_Of_Frame_ORFs))
        per_Expanded_CDS = comp.expanded_CDS * 100 / len(comp.genes_Detected)
        per_Expanded_Start = comp.expanded_Start * 100 / len(comp.genes_Detected)
        per_Expanded_Stop = comp.expanded_Stop * 100 / len(comp.genes_Detected)
        perfect_Matches_Percentage = comp.perfect_Matches * 100 / len(comp.matched_ORFs)
        perfect_Starts_Percentage = float(comp.perfect_Starts) * float(100) / float(len(comp.matched_ORFs))
        perfect_Stops_Percentage = float(comp.perfect_Stops) * float(100) / float(len(comp.matched_ORFs))
    except ZeroDivisionError:
        correct_Frame_Percentage = 0
        per_Expanded_CDS = 0
        per_Expanded_Start = 0
        per_Expanded_Stop = 0
        perfect_Matches_Percentage = 0
        perfect_Starts_Percentage = 0
        perfect_Stops_Percentage = 0
    #############################################
    #Missed Genes  Metrics:
    mg_Starts = []
    mg_Stops = []
    mg_Lengths = []
    mg_Strands = []
    for mg, seq in comp.genes_Undetected.items():
        mg = mg.split(',')
        mg_Starts.append(mg[3])
        mg_Stops.append(mg[4])
        mg_Strands.append(mg[2])
        mg_Lengths.append(int(mg[1])-int(mg[0]))

    mg_ATG = mg_Starts.count('ATG') *100 / len(comp.genes_Undetected)
    mg_GTG = mg_Starts.count('GTG') *100 / len(comp.genes_Undetected)
    mg_TTG = mg_Starts.count('TTG') *100 / len(comp.genes_Undetected)
    mg_ATT = mg_Starts.count('ATT') *100 / len(comp.genes_Undetected)
    mg_CTG = mg_Starts.count('CTG') *100 / len(comp.genes_Undetected)
    mg_O_Start = 100 - (mg_ATG+mg_GTG+mg_TTG+mg_ATT+mg_CTG)
    mg_TGA = mg_Stops.count('TGA') *100 / len(comp.genes_Undetected)
    mg_TAA = mg_Stops.count('TAA') *100 / len(comp.genes_Undetected)
    mg_TAG = mg_Stops.count('TAG') *100 / len(comp.genes_Undetected)
    mg_O_Stop = 100 - (mg_TGA+mg_TAA+mg_TAG)
    median_mg_Len = np.median(mg_Lengths)
    mg_Pos = mg_Strands.count('+')
    mg_Neg = mg_Strands.count('-')
    Missed_Gene_Metrics = (format(mg_ATG,'.2f'),format(mg_GTG,'.2f'),format(mg_TTG,'.2f'),format(mg_ATT,'.2f'),format(mg_CTG,'.2f'),format(mg_O_Start,'.2f'),format(mg_TGA,'.2f'),format(mg_TAA,'.2f'),format(mg_TAG,'.2f'),format(mg_O_Stop,'.2f'),format(median_mg_Len,'.2f'),mg_Pos,mg_Neg)

    # Unmathced ORF Metrics:
    uo_Starts = []
    uo_Stops = []
    uo_Lengths = []
    uo_Strands = []
    for uo, seq in comp.unmatched_ORFs.items():
        uo = uo.split(',')
        uo_Starts.append(uo[3])
        uo_Stops.append(uo[4])
        uo_Strands.append(uo[2])
        uo_Lengths.append(int(uo[1]) - int(uo[0]))
    uo_ATG = uo_Starts.count('ATG') * 100 / len(comp.unmatched_ORFs)
    uo_GTG = uo_Starts.count('GTG') * 100 / len(comp.unmatched_ORFs)
    uo_TTG = uo_Starts.count('TTG') * 100 / len(comp.unmatched_ORFs)
    uo_ATT = uo_Starts.count('ATT') * 100 / len(comp.unmatched_ORFs)
    uo_CTG = uo_Starts.count('CTG') * 100 / len(comp.unmatched_ORFs)
    uo_O_Start = 100 - (uo_ATG + uo_GTG + uo_TTG + uo_ATT + uo_CTG)
    uo_TGA = uo_Stops.count('TGA') * 100 / len(comp.unmatched_ORFs)
    uo_TAA = uo_Stops.count('TAA') * 100 / len(comp.unmatched_ORFs)
    uo_TAG = uo_Stops.count('TAG') * 100 / len(comp.unmatched_ORFs)
    uo_O_Stop = len(uo_Stops) -(uo_TGA+uo_TAA+uo_TAG)
    uo_O_Stop = uo_O_Stop *100 / len(comp.unmatched_ORFs)
    median_uo_Len = np.median(uo_Lengths)
    uo_Pos = uo_Strands.count('+')
    uo_Neg = uo_Strands.count('-')
    unmatched_orf_metrics = (format(uo_ATG,'.2f'),format(uo_GTG,'.2f'),format(uo_TTG,'.2f'),format(uo_ATT,'.2f'),format(uo_CTG,'.2f'),format(uo_O_Start,'.2f'),format(uo_TGA,'.2f'),format(uo_TAA,'.2f'),format(uo_TAG,'.2f'),format(uo_O_Stop,'.2f'),format(median_uo_Len,'.2f'),uo_Pos,uo_Neg)
    #################################
    metric_description = ['Number of ORFs',	'Percentage Difference of ORFs', 'Number of ORFs that Detected a Gene', 'Percentage of ORFs that Detected a Gene', 'Number of Genes Detected',
                        'Percentage of Genes Detected', 'Median Length of All ORFs', 'Median Length Difference', 'Minimum Length of All ORFs', 'Minimum Length Difference',
                        'Maximum Length of All ORFs', 'Maximum Length Difference', 'Number of Perfect Matches', 'Percentage of Perfect Matches', 'Number of Perfect Starts',
                        'Percentage of Perfect Starts', 'Number of Perfect Stops',	'Percentage of Perfect Stops', 'Number of Matched ORFs in Correct Frame', 'Percentage of Matched ORFs in Correct Frame',
                        'Number of Matched ORFs Expanding a Coding Region', 'Percentage of Matched ORFs Expanding a Coding Region', 'Percentage of Matched ORFs Expanding Start Region',
                        'Percentage of Matched ORFs Expanding Stop Region', 'Number of All ORFs on Positive Strand', 'Percentage of All ORFs on Positive Strand',
                        'Number of All ORFs on Negative Strand', 'Percentage of All ORFs on Negative Strand', 'Median Start Difference of Matched ORFs', 'Median Stop Difference of Matched ORFs','ATG Start Percentage',
                        'GTG Start Percentage', 'TTG Start Percentage', 'ATT Start Percentage', 'CTG Start Percentage', 'Other Start Codon Percentage', 'TAG Stop Percentage', 'TAA Stop Percentage',
                        'TGA Stop Percentage', 'Other Stop Codon Percentage', 'True Positive', 'False Positive', 'False Negative', 'Precision', 'Recall', 'False Discovery Rate',
                        'Nucleotide True Positive', 'Nucleotide False Positive', 'Nucleotide True Negative', 'Nucleotide False Negative', 'Nucleotide Precision', 'Nucleotide Recall',
                        'Nucleotide False Discovery Rate','ORF Nucleotide Coverage of Genome','Matched ORF Nucleotide Coverage of Genome']

    metrics = [len(orfs), format(ORFs_Diff,'.2f'), len(comp.matched_ORFs), format(matched_ORF_Percentage,'.2f'), len(comp.genes_Detected),
              format(genes_Detected_Percentage,'.2f'), format(median_ORF_Length,'.2f'), format(median_Length_Diff,'.2f'), min_ORF_Length, format(min_Length_Diff,'.2f'),
              max_ORF_Length, format(max_Length_Diff,'.2f'),comp.perfect_Matches, format(perfect_Matches_Percentage,'.2f'),  comp.perfect_Starts,
              format(perfect_Starts_Percentage,'.2f'), comp.perfect_Stops, format(perfect_Stops_Percentage,'.2f'), comp.correct_Frame_Number, format(correct_Frame_Percentage,'.2f'),
              format(per_Expanded_CDS,'.2f'), format(per_Expanded_Start,'.2f'), format(per_Expanded_Stop,'.2f'), comp.pos_Strand, format(pos_Strand_Percentage,'.2f'), comp.neg_Strand, format(neg_Strand_Percentage,'.2f'),
              format(median_Start_Difference,'.2f'), format(median_Stop_Difference,'.2f'), format(atg_P, '.2f'), format(gtg_P,'.2f'), format(ttg_P, '.2f'), format(att_P, '.2f'), format(ctg_P, '.2f'), format(other_Start_P, '.2f'),
              format(tag_P, '.2f'), format(taa_P, '.2f'), format(tga_P, '.2f'), format(other_Stop_P,'.2f'), format(TP,'.2f'), format(FP,'.2f'), format(FN,'.2f'), format(precision,'.2f'), format(recall,'.2f'),
              format(false_Discovery_Rate,'.2f'), format(NT_TP,'.2f'), format(NT_FP,'.2f'), format(NT_TN,'.2f'), format(NT_FN,'.2f'), format(NT_Precision,'.2f'), format(NT_Recall,'.2f'),
              format(NT_False_Discovery_Rate,'.2f'),format(orf_Coverage_Genome,'.2f'),format(matched_ORF_Coverage_Genome,'.2f')]

    rep_metric_description = ['Percentage of ORFs that Detected a Gene','Percentage Difference of ORFs','Median Length Difference','Percentage of Perfect Matches','Percentage of Perfect Starts','Percentage of Perfect Stops',
                               'Median Start Difference of Matched ORFs', 'Median Stop Difference of Matched ORFs','Precision', 'Recall']

    rep_metrics = [format(matched_ORF_Percentage,'.2f'),format(ORFs_Diff,'.2f'),format(median_Length_Diff,'.2f'),format(perfect_Matches_Percentage,'.2f'), format(perfect_Starts_Percentage,'.2f'), format(perfect_Stops_Percentage,'.2f'),
                   format(median_Start_Difference,'.2f'), format(median_Stop_Difference,'.2f'),format(precision,'.2f'), format(recall,'.2f')]

    return metric_description, metrics, rep_metric_description, rep_metrics, start_Difference, stop_Difference, other_Starts, other_Stops, comp.genes_Undetected, comp.unmatched_ORFs, Missed_Gene_Metrics, unmatched_orf_metrics,gene_Coverage_Genome
