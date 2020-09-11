import collections

import numpy as np
from Tools.utils import revCompIterative


class comparator:  # Class to hold global-type variables
    def __init__(self, perfect_Starts=0, perfect_Stops=0, perfect_Matches=0, genome_Seq='',
                 genome_Seq_Rev='',
                 genome_Size=0, correct_Frame_Number=0, extended_Start=0,
                 extended_Stop=0, extended_CDS=0, matched_ORFs=collections.OrderedDict(),multi_Matched_ORFs=collections.OrderedDict(),
                 unmatched_ORFs=collections.OrderedDict(), genes_Detected=collections.OrderedDict(),
                 genes_Undetected=collections.OrderedDict(),
                 out_Of_Frame_ORFs=collections.OrderedDict(), start_Difference=[], stop_Difference=[],
                 orf_Lengths=[], gene_Lengths=[], gene_Pos_Olap=[], gene_Neg_Olap=[], orf_Pos_Olap=[], orf_Neg_Olap=[], m_ORF_Pos_Olap=[], m_ORF_Neg_Olap=[], gene_GC=[],
                 orf_GC=[], m_ORF_GC=[], gene_Short=[], orf_Short=[], m_ORF_Short = [], pos_Strand=0, neg_Strand=0):
        self.perfect_Starts, self.perfect_Stops, self.perfect_Matches, self.genome_Seq, self.genome_Seq_Rev, self.genome_Size, self.correct_Frame_Number, self.extended_Start, self.extended_Stop, self.extended_CDS, \
        self.matched_ORFs, self.multi_Matched_ORFs, self.unmatched_ORFs, self.genes_Detected, self.genes_Undetected, self.out_Of_Frame_ORFs, self.start_Difference, self.stop_Difference, self.orf_Lengths, self.gene_Lengths, self.gene_Pos_Olap, \
        self.gene_Neg_Olap, self.orf_Pos_Olap, self.orf_Neg_Olap, self.m_ORF_Pos_Olap, self.m_ORF_Neg_Olap, self.gene_GC, self.orf_GC, self.m_ORF_GC, self.gene_Short, self.orf_Short, self.m_ORF_Short, self.pos_Strand, \
        self.neg_Strand = perfect_Starts, perfect_Stops, perfect_Matches, genome_Seq, genome_Seq_Rev, \
        genome_Size, correct_Frame_Number, extended_Start, extended_Stop, extended_CDS, matched_ORFs, multi_Matched_ORFs, unmatched_ORFs, genes_Detected, genes_Undetected, out_Of_Frame_ORFs, start_Difference, stop_Difference, orf_Lengths, \
        gene_Lengths, gene_Pos_Olap, gene_Neg_Olap, orf_Pos_Olap, orf_Neg_Olap, m_ORF_Pos_Olap, m_ORF_Neg_Olap, gene_GC, orf_GC, m_ORF_GC,  gene_Short, orf_Short, m_ORF_Short, pos_Strand, neg_Strand

comp = comparator()





def nuc_Count(start,stop,strand): # Gets correct seq then returns GC
    if strand == '-':
        r_Start = comp.genome_Size - stop
        r_Stop = comp.genome_Size - start
        seq = (comp.genome_Seq_Rev[r_Start:r_Stop + 1])
    elif strand == '+':
        seq = (comp.genome_Seq[start - 1:stop])
    c = 0
    a = 0
    g = 0
    t = 0
    n = 0
    for i in seq:
        if "C" in i:
            c += 1
        elif "G" in i:
            g += 1
        elif "A" in i:
            a += 1
        elif "T" in i:
            t+=1
        elif "N" in i:
            n+=1
    gc_content = (g + c) * 100 / (a + t + g + c + n)
    n_per = n * 100 / (a + t + g + c + n)
    return gc_content


def orf_Unmatched(o_Start, o_Stop, o_Strand):
    if o_Strand == '-':
        r_Start = comp.genome_Size - o_Stop
        r_Stop = comp.genome_Size - o_Start
        Unmatched_ORF = str(o_Start) + ',' + str(o_Stop) + ',' + o_Strand + ',' +comp.genome_Seq_Rev[r_Start:r_Start+3] + ',' +comp.genome_Seq_Rev[r_Stop - 2:r_Stop + 1]
        seq = (comp.genome_Seq_Rev[r_Start:r_Stop+1])
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
        comp.extended_CDS +=1
    if o_Start < g_Start:
        comp.extended_Start +=1
    if o_Stop > g_Stop:
        comp.extended_Stop +=1

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
    atg_P = format(100* atg / len(orfs),'.2f')
    gtg_P = format(100 * gtg / len(orfs),'.2f')
    ttg_P = format(100 * ttg / len(orfs),'.2f')
    att_P = format(100 * att  / len(orfs),'.2f')
    ctg_P = format(100 * ctg  / len(orfs),'.2f')
    other_Start_P = format(100 * other / len(orfs),'.2f')
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
    tag_p = format(100 * tag  / len(orfs),'.2f')
    taa_p = format(100 * taa  / len(orfs),'.2f')
    tga_p = format(100 * tga / len(orfs),'.2f')
    other_Stop_P = format(100 * other / len(orfs),'.2f')
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
    for gene_Num, gene_Details in genes.items(): #Loop through each gene to compare against predicted ORFs
        gene_Details = gene_Details.split(',')
        g_Start = int(gene_Details[0])
        g_Stop = int(gene_Details[1])
        g_Strand = gene_Details[2]
        g_pos = str(g_Start)+','+str(g_Stop)
        gene_Set = set(range(g_Start, g_Stop + 1)) # Used to check Overlap of ORFs and pick best match - slow but confirms best match
        overlapping_ORFs = collections.OrderedDict()
        perfect_Match = False
        out_Frame = False
        for pos, orf_Details  in orfs.items(): #Check if perfect match, if not check if match covers at least 75% of gene - Loop through ALL ORFs - SLOW
            o_Start = int(pos.split(',')[0])
            o_Stop = int(pos.split(',')[1])
            o_Strand = orf_Details[0]
            orf_Set = set(range(o_Start, o_Stop + 1))
            if o_Stop <= g_Start or o_Start >= g_Stop: #Not caught up yet or ORFs may not be ordered - Slow but needed for unordered ORF Predictions
                continue
            elif o_Start == g_Start and o_Stop == g_Stop: #If perfect match, break and skip the rest of the ORFs
                perfect_Match = True
                break
            elif g_Start <= o_Start < g_Stop or g_Start < o_Stop < g_Stop: # If ORF Start or Stop is between gene Start or Stop
                overlap = len(gene_Set.intersection(orf_Set))
                coverage = 100 * float(overlap) / float(len(gene_Set))
                orf_Details.append(coverage)
                if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and coverage >= 75 :  #Only continue if ORF covers at least 75% of the gene and is in frame
                    overlapping_ORFs.update({pos:orf_Details})
                elif coverage >= 75: # Not in frame / on same strand
                    comp.out_Of_Frame_ORFs.update({pos:orf_Details})
                    out_Frame = True
            elif o_Start <= g_Start and o_Stop >= g_Stop: # If ORF extends one or both ends of the gene
                overlap = len(gene_Set.intersection(orf_Set))
                coverage = 100 * float(overlap) / float(len(gene_Set))
                orf_Details.append(coverage)
                if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and coverage >= 75:  #Only continue if ORF covers at least 75% of the gene and is in frame
                    overlapping_ORFs.update({pos:orf_Details})
                elif coverage >= 75:
                    comp.out_Of_Frame_ORFs.update({pos:orf_Details})
                    out_Frame = True
            else:
                print("Unexpected Error Finding ORFs") # Should not happen
        #Now Check that we select the best ORF
                                                                                                                        ### Multi_Match_ORFs Should contain All genes found by a specific ORF
        if perfect_Match == True: # Check if the ORF is a perfect match to the Gene
            m_ORF_Details = orf_Details
            m_ORF_Details.append(g_pos)
            if g_pos in comp.matched_ORFs.keys():
                comp.multi_Matched_ORFs.update(({g_pos:m_ORF_Details}))
            comp.matched_ORFs.update({g_pos:m_ORF_Details})
            comp.genes_Detected.update({str(gene_Details):g_pos})
            orf_Details.append(100)
            match_Statistics(o_Start, o_Stop, g_Start, g_Stop)
            comp.perfect_Matches += 1
            print('Perfect Match')
        elif perfect_Match == False and len(overlapping_ORFs) == 1: # If we do not have a perfect match but 1 ORF which has passed the filtering
            orf_Pos = list(overlapping_ORFs.keys())[0]
            o_Start = int(orf_Pos.split(',')[0])
            o_Stop = int(orf_Pos.split(',')[1])
            orf_Details = overlapping_ORFs[orf_Pos]
            m_ORF_Details = orf_Details
            m_ORF_Details.append(g_pos)
            if orf_Pos in comp.matched_ORFs.keys():
                comp.multi_Matched_ORFs.update(({g_pos:m_ORF_Details}))
            comp.matched_ORFs.update({orf_Pos:m_ORF_Details})
            comp.genes_Detected.update({str(gene_Details):orf_Pos})
            match_Statistics(o_Start,o_Stop,g_Start,g_Stop)
            print('Partial Match')
        elif perfect_Match == False and len(overlapping_ORFs) >= 1: # If we have more than 1 potential ORF match, we check to see which is the 'best' hit
            orf_Pos,orf_Details = candidate_ORF_Selection(gene_Set,overlapping_ORFs) # Return best match
            o_Start = int(orf_Pos.split(',')[0])
            o_Stop = int(orf_Pos.split(',')[1])
            m_ORF_Details = orf_Details
            m_ORF_Details.append(g_pos)
            if orf_Pos in comp.matched_ORFs.keys():
                comp.multi_Matched_ORFs.update(({g_pos:m_ORF_Details}))
            comp.matched_ORFs.update({orf_Pos:m_ORF_Details})
            comp.genes_Detected.update({str(gene_Details):orf_Pos})
            match_Statistics(o_Start,o_Stop,g_Start,g_Stop)
            print('There was more than 1 potential Match - Best Chosen')
        elif out_Frame: # Keep record of ORFs which overlap a gene but in the wrong frame
            print("Out of Frame ORF")
            genes_Unmatched(g_Start, g_Stop, g_Strand) #
        else:
            genes_Unmatched(g_Start, g_Stop, g_Strand) # No hit
            print("No Hit")
    for key in comp.matched_ORFs: # Remove ORFs from out of frame if ORF was correctly matched to another Gene
        if key in comp.out_Of_Frame_ORFs:
            del comp.out_Of_Frame_ORFs[key]

    ######################################## ORF Lengths and Precision
    start_Difference = [x for x in comp.start_Difference if x != 0] # Remove 0s (Perfect hits)
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
    gene_Nuc_Array = np.zeros((comp.genome_Size), dtype=np.int)
    orf_Nuc_Array = np.zeros((comp.genome_Size), dtype=np.int)
    matched_ORF_Nuc_Array = np.zeros((comp.genome_Size), dtype=np.int)

    gene_Prev_Stop = 0
    for gene_Num, gene_Details in genes.items(): #Loop through each gene to compare against predicted ORFs
        gene_Details = gene_Details.split(',')
        g_Start = int(gene_Details[0])
        g_Stop = int(gene_Details[1])
        g_Strand = gene_Details[2]
        gene_Length = (g_Stop - g_Start)
        comp.gene_Lengths.append(gene_Length)
        gene_Nuc_Array[g_Start-1:g_Stop] = [1] # Changing all between the two positions to 1's
        if gene_Prev_Stop > g_Start: #Check if prev gene overlaps current gene
            if '+' in g_Strand:
                comp.gene_Pos_Olap.append(gene_Prev_Stop - g_Start)
            elif '-' in g_Strand:
                comp.gene_Neg_Olap.append(gene_Prev_Stop - g_Start)
        gene_Prev_Stop = g_Stop
        if gene_Length < 100:
            comp.gene_Short.append(gene_Length)
        comp.gene_GC.append(nuc_Count(g_Start, g_Stop, g_Strand))
    min_Gene_Length = min(comp.gene_Lengths)
    max_Gene_Length = max(comp.gene_Lengths)
    median_Gene_Length = np.median(comp.gene_Lengths)

    orf_Prev_Stop = 0
    for o_Positions,orf_Details in orfs.items():
        o_Start = int(o_Positions.split(',')[0])
        o_Stop = int(o_Positions.split(',')[1])
        o_Strand = orf_Details[0]
        orf_Length = (o_Stop - o_Start)
        comp.orf_Lengths.append(orf_Length)
        orf_Nuc_Array[o_Start-1:o_Stop] = [1] # Changing all between the two positions to 1's
        if orf_Prev_Stop > o_Start: #Check if prev orf overlaps current orf
            if '+' in o_Strand:
                comp.orf_Pos_Olap.append(orf_Prev_Stop - o_Start)
            elif '-' in o_Strand:
                comp.orf_Neg_Olap.append(orf_Prev_Stop - o_Start)
        orf_Prev_Stop = o_Stop
        if orf_Length < 100:
            comp.orf_Short.append(orf_Length)
        comp.orf_GC.append(nuc_Count(o_Start,o_Stop,o_Strand))

        # Get ORF Strand metrics:
        if o_Strand == "+":  # Get number of Positive and Negative strand ORFs
            comp.pos_Strand += 1
        elif o_Strand == "-":
            comp.neg_Strand += 1
        # Stats just for Unmatched ORFs
        if o_Positions not in list(comp.matched_ORFs.keys()):
            orf_Unmatched(o_Start, o_Stop, o_Strand)
    # Nucleotide Coverage calculated from ORFs matching a gene only
    matched_ORF_Prev_Stop = 0
    for mo_Positions, m_ORF_Details in comp.matched_ORFs.items():
        mo_Start = int(mo_Positions.split(',')[0])
        mo_Stop = int(mo_Positions.split(',')[1])
        mo_Strand = m_ORF_Details[0]
        mo_Length = (mo_Stop - mo_Start)
        matched_ORF_Nuc_Array[mo_Start-1 :mo_Stop] = [1] # Changing all between the two positions to 1's
        if matched_ORF_Prev_Stop > mo_Start: #Check if prev orf overlaps current orf
            if '+' in mo_Strand:
                comp.m_ORF_Pos_Olap.append(matched_ORF_Prev_Stop - mo_Start)
            elif '-' in mo_Strand:
                comp.m_ORF_Neg_Olap.append(matched_ORF_Prev_Stop - mo_Start)
        matched_ORF_Prev_Stop = mo_Stop
        if mo_Length < 100:
            comp.m_ORF_Short.append(mo_Length)
        comp.m_ORF_GC.append(nuc_Count(mo_Start,mo_Stop,mo_Strand))


    gene_Coverage_Genome = format(100 * np.count_nonzero(gene_Nuc_Array) / comp.genome_Size,'.2f')
    orf_Coverage_Genome = format(100 * np.count_nonzero(orf_Nuc_Array) / comp.genome_Size,'.2f')
    matched_ORF_Coverage_Genome = format(100 * np.count_nonzero(matched_ORF_Nuc_Array) / comp.genome_Size,'.2f')
    # gene and orf nucleotide Intersection
    gene_count = np.count_nonzero(gene_Nuc_Array)
    orf_count = np.count_nonzero(orf_Nuc_Array)
    gene_ORF_Nuc_Intersection = np.count_nonzero(gene_Nuc_Array & orf_Nuc_Array)
    #not gene but orf nucleotides
    not_Gene_Nuc_Array = np.logical_not(gene_Nuc_Array) + [0 for i in range(len(gene_Nuc_Array))] #End part to keep array as 1,0 not T,F
    not_Gene_Nuc_And_ORF_Count = np.count_nonzero(not_Gene_Nuc_Array & orf_Nuc_Array)
    #not orf nucleotides but gene
    not_ORF_Nuc_Array = np.logical_not(orf_Nuc_Array) + [0 for i in range(len(orf_Nuc_Array))] #End part to keep array as 1,0 not T,F
    not_ORF_Nuc_And_Gene_Count = np.count_nonzero(not_ORF_Nuc_Array & gene_Nuc_Array)
    #not gene or orf nucleotides
    not_Gene_Nuc_Not_ORF_Nuc_Count = np.count_nonzero(not_Gene_Nuc_Array & not_ORF_Nuc_Array)
    #Nucleotide 'accuracy' - Normalised by number of nucelotides annotated by a gene
    NT_TP = format(gene_ORF_Nuc_Intersection  / np.count_nonzero(gene_Nuc_Array),'.2f')
    NT_FP = format(not_Gene_Nuc_And_ORF_Count / np.count_nonzero(not_Gene_Nuc_Array),'.2f')
    NT_FN = format(not_ORF_Nuc_And_Gene_Count / np.count_nonzero(gene_Nuc_Array),'.2f')
    NT_TN = format(not_Gene_Nuc_Not_ORF_Nuc_Count / np.count_nonzero(not_Gene_Nuc_Array),'.2f')
    NT_Precision = format(gene_ORF_Nuc_Intersection / (gene_ORF_Nuc_Intersection + not_Gene_Nuc_And_ORF_Count),'.2f')
    NT_Recall = format(gene_ORF_Nuc_Intersection / (gene_ORF_Nuc_Intersection + not_ORF_Nuc_And_Gene_Count),'.2f')
    NT_False_Discovery_Rate = format(not_Gene_Nuc_And_ORF_Count / (not_Gene_Nuc_And_ORF_Count + gene_ORF_Nuc_Intersection),'.2f')
    ################################# Precision and Recall of whole ORFs and Genes
    TP = format(len(comp.genes_Detected)  / len(genes),'.2f')
    FP = format(len(comp.unmatched_ORFs) / len(genes),'.2f')
    FN = format(len(comp.genes_Undetected)  / len(genes),'.2f')
    #try: # Incase no ORFs found a gene
    precision = format(float(TP)/(float(TP)+float(FP)),'.2f')
    recall = format(float(TP)/(float(TP)+float(FN)),'.2f')
    false_Discovery_Rate = format(float(FP)/(float(FP)+float(TP)),'.2f')
    min_ORF_Length = min(comp.orf_Lengths)
    max_ORF_Length = max(comp.orf_Lengths)
    median_ORF_Length = np.median(comp.orf_Lengths)
    # except ZeroDivisionError:
    #     precision = 0
    #     recall = 0
    #     false_Discovery_Rate = 0
    #     min_ORF_Length = 0
    #     max_ORF_Length = 0
    #     median_ORF_Length = 0


##########################################################################
    # Metrics - There are numerous cases where certain metric calculations may return a ZeroDivError.
    ORFs_Difference = format(100 * (len(orfs) - len(genes)) / len(genes),'.2f') #Difference off +/-
    genes_Detected_Percentage = format(100 * (len(comp.genes_Detected) / len(genes)),'.2f')
    matched_ORF_Percentage = format(100 * (len(comp.matched_ORFs) / len(orfs)),'.2f')
    all_ORF_Olap = (comp.orf_Pos_Olap + comp.orf_Neg_Olap) #Combine pos and neg strand overlaps
    matched_ORF_Olap = (comp.m_ORF_Pos_Olap + comp.m_ORF_Neg_Olap)
    all_Gene_Olap = (comp.gene_Pos_Olap + comp.gene_Neg_Olap)

    if all_ORF_Olap: # If no overlapping ORFs
        overlap_Difference = format(100 * (len(all_ORF_Olap) - len(all_Gene_Olap)) / len(all_Gene_Olap),'.2f')
        matched_Overlap_Difference = format(100 * (len(matched_ORF_Olap) - len(all_Gene_Olap)) / len(all_Gene_Olap),'.2f')
        num_All_ORF_Olap = len(all_ORF_Olap)
        if matched_ORF_Olap:
            max_Matched_ORF_Olap = max(matched_ORF_Olap)
            matched_Median_ORF_Overlap = format(np.median(matched_ORF_Olap), '.2f')
        else:
            max_Matched_ORF_Olap = 'N/A'
            matched_Median_ORF_Overlap = 'N/A'
        max_All_ORF_Olap = max(all_ORF_Olap)
        median_ORF_Overlap = format(np.median(all_ORF_Olap), '.2f')

    else:
        overlap_Difference = -100
        matched_Overlap_Difference = 'N/A'
        num_All_ORF_Olap = 0
        max_Matched_ORF_Olap = 'N/A'
        max_All_ORF_Olap = 'N/A'
        median_ORF_Overlap = 'N/A'
        matched_Median_ORF_Overlap = 'N/A'


    if comp.orf_Short and comp.gene_Short: # IF Short-ORFs/Genes
        short_ORF_Difference = format(100 * (len(comp.orf_Short) - len(comp.gene_Short))/ len(comp.gene_Short),'.2f')
        matched_Short_ORF_Difference = format(100 * (len(comp.m_ORF_Short) - len(comp.gene_Short))/ len(comp.gene_Short),'.2f')
        num_ORF_Short = len(comp.orf_Short)
        num_Matched_ORF_Short = len(comp.m_ORF_Short)
    elif comp.orf_Short: # If only Short-ORFs
        num_ORF_Short = len(comp.orf_Short)
        num_Matched_ORF_Short = len(comp.m_ORF_Short)
        short_ORF_Difference = (num_ORF_Short*100)
        matched_Short_ORF_Difference = 'N/A'
    elif comp.gene_Short: # If only Short-Genes
        short_ORF_Difference = -100
        matched_Short_ORF_Difference = 'N/A'
        num_ORF_Short = 0
        num_Matched_ORF_Short = 'N/A'


    median_Length_Difference = format(100 * (median_ORF_Length - median_Gene_Length) / median_Gene_Length,'.2f')
    min_Length_Difference = format(100 * (min_ORF_Length - min_Gene_Length) / min_Gene_Length,'.2f')
    max_Length_Difference = format(100 * (max_ORF_Length - max_Gene_Length) / max_Gene_Length,'.2f')
    pos_Strand_Percentage = format(comp.pos_Strand / len(orfs),'.2f')
    neg_Strand_Percentage = format(comp.neg_Strand / len(orfs),'.2f')
    median_ORF_GC = format(np.median(comp.orf_GC),'.2f')
    matched_Median_ORF_GC = format(np.median(comp.m_ORF_GC),'.2f')
    median_Gene_GC = format(np.median(comp.gene_GC),'.2f')
    median_GC_Difference = format(100 * (float(median_ORF_GC) - float(median_Gene_GC)) / float(median_Gene_GC),'.2f')
    matched_Median_GC_Difference = format(100 * (float(matched_Median_ORF_GC) - float(median_Gene_GC)) / float(median_Gene_GC),'.2f')
    #############################################
    if comp.matched_ORFs: # Incase no ORFs detected a gene
        #correct_Frame_Number = comp.correct_Frame_Number * 100 / (len(comp.genes_Detected) + len(comp.out_Of_Frame_ORFs))
        extended_CDS_Percentage = format(100 * comp.extended_CDS / len(comp.matched_ORFs),'.2f')
        extended_Start_Percentage = format(100 * comp.extended_Start / len(comp.matched_ORFs),'.2f')
        extended_Stop_Percentage = format(100 * comp.extended_Stop / len(comp.matched_ORFs),'.2f')
        perfect_Matches_Percentage = format(100 * comp.perfect_Matches / len(comp.matched_ORFs),'.2f')
        perfect_Starts_Percentage = format(100 * comp.perfect_Starts / len(comp.matched_ORFs),'.2f')
        perfect_Stops_Percentage = format(100 * comp.perfect_Stops / len(comp.matched_ORFs),'.2f')
    else:
        #correct_Frame_Percentage = 0
        extended_CDS_Percentage = 0
        extended_Start_Percentage = 0
        extended_Stop_Percentage = 0
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

    mg_ATG = 100 * mg_Starts.count('ATG')  / len(comp.genes_Undetected)
    mg_GTG = 100 * mg_Starts.count('GTG')  / len(comp.genes_Undetected)
    mg_TTG = 100 * mg_Starts.count('TTG')  / len(comp.genes_Undetected)
    mg_ATT = 100 * mg_Starts.count('ATT')  / len(comp.genes_Undetected)
    mg_CTG = 100 * mg_Starts.count('CTG')  / len(comp.genes_Undetected)
    mg_O_Start = 100 - (mg_ATG+mg_GTG+mg_TTG+mg_ATT+mg_CTG)
    mg_TGA = 100 * mg_Stops.count('TGA')  / len(comp.genes_Undetected)
    mg_TAA = 100 * mg_Stops.count('TAA')  / len(comp.genes_Undetected)
    mg_TAG = 100 * mg_Stops.count('TAG')  / len(comp.genes_Undetected)
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
    uo_ATG = 100 * uo_Starts.count('ATG') / len(comp.unmatched_ORFs)
    uo_GTG = 100 * uo_Starts.count('GTG') / len(comp.unmatched_ORFs)
    uo_TTG = 100 * uo_Starts.count('TTG') / len(comp.unmatched_ORFs)
    uo_ATT = 100 * uo_Starts.count('ATT') / len(comp.unmatched_ORFs)
    uo_CTG = 100 * uo_Starts.count('CTG') / len(comp.unmatched_ORFs)
    uo_O_Start = 100 - (uo_ATG + uo_GTG + uo_TTG + uo_ATT + uo_CTG)
    uo_TGA = 100 * uo_Stops.count('TGA') / len(comp.unmatched_ORFs)
    uo_TAA = 100 * uo_Stops.count('TAA') / len(comp.unmatched_ORFs)
    uo_TAG = 100 * uo_Stops.count('TAG') / len(comp.unmatched_ORFs)
    uo_O_Stop = len(uo_Stops) -(uo_TGA+uo_TAA+uo_TAG)
    uo_O_Stop = 100 * uo_O_Stop  / len(comp.unmatched_ORFs)
    median_uo_Len = np.median(uo_Lengths)
    uo_Pos = uo_Strands.count('+')
    uo_Neg = uo_Strands.count('-')
    unmatched_orf_metrics = (format(uo_ATG,'.2f'),format(uo_GTG,'.2f'),format(uo_TTG,'.2f'),format(uo_ATT,'.2f'),format(uo_CTG,'.2f'),format(uo_O_Start,'.2f'),format(uo_TGA,'.2f'),format(uo_TAA,'.2f'),format(uo_TAG,'.2f'),format(uo_O_Stop,'.2f'),format(median_uo_Len,'.2f'),uo_Pos,uo_Neg)
    #################################



    all_Metrics = collections.OrderedDict({'Number of ORFs':len(orfs),'Percent Difference of All ORFs':ORFs_Difference,'Number of ORFs that Detected a Gene':len(comp.matched_ORFs),'Percentage of ORFs that Detected a Gene':matched_ORF_Percentage,
                                           'Number of Genes Detected':len(comp.genes_Detected),'Percentage of Genes Detected':genes_Detected_Percentage,'Median Length of All ORFs':median_ORF_Length,'Median Length Difference':median_Length_Difference,
                                           'Minimum Length of All ORFs':min_ORF_Length,'Minimum Length Difference':min_Length_Difference,'Maximum Length of All ORFs':max_ORF_Length,'Maximum Length Difference':max_Length_Difference,
                                           'Median GC content of All ORFs':median_ORF_GC,'Percent Difference of All ORFs Median GC':median_GC_Difference,'Median GC content of Matched ORFs':matched_Median_ORF_GC,
                                           'Percent Difference of Matched ORF GC':matched_Median_GC_Difference,'Number of ORFs which Overlap Another ORF':num_All_ORF_Olap,'Percent Difference of Overlapping ORFs':overlap_Difference,
                                           'Maximum ORF Overlap':max_All_ORF_Olap,'Median ORF Overlap':median_ORF_Overlap,'Number of Matched ORFs Overlapping Another ORF':len(matched_ORF_Olap),'Percent Difference of Matched ORFs Overlapping Another ORF':matched_Overlap_Difference,
                                           'Maximum Matched ORF Overlap':max_Matched_ORF_Olap,'Median Matched ORF Overlap':matched_Median_ORF_Overlap,'Number of Short-ORFs':num_ORF_Short,'Percent Difference of Short-ORFs':short_ORF_Difference,
                                           'Number of Short-Matched-ORFs':num_Matched_ORF_Short,'Percent Difference of Short-Matched-ORFs':matched_Short_ORF_Difference,'Number of Perfect Matches':comp.perfect_Matches,'Percentage of Perfect Matches':perfect_Matches_Percentage,
                                           'Number of Perfect Starts':comp.perfect_Starts,'Percentage of Perfect Starts':perfect_Starts_Percentage,'Number of Perfect Stops':comp.perfect_Stops,'Percentage of Perfect Stops':perfect_Stops_Percentage,
                                           'Number of Out of Frame ORFs':len(comp.out_Of_Frame_ORFs),'Number of Matched ORFs Extending a Coding Region':comp.extended_CDS,'Percentage of Matched ORFs Extending a Coding Region':extended_CDS_Percentage,
                                           'Number of Matched ORFs Extending Start Region':comp.extended_Start,'Percentage of Matched ORFs Extending Start Region':extended_Start_Percentage,'Number of Matched ORFs Extending Stop Region':comp.extended_Stop,
                                           'Percentage of Matched ORFs Extending Stop Region':extended_Stop_Percentage,'Number of All ORFs on Positive Strand':comp.pos_Strand,'Percentage of All ORFs on Positive Strand':pos_Strand_Percentage,
                                           'Number of All ORFs on Negative Strand':comp.neg_Strand,'Percentage of All ORFs on Negative Strand':neg_Strand_Percentage,'Median Start Difference of Matched ORFs':median_Start_Difference,
                                           'Median Stop Difference of Matched ORFs':median_Stop_Difference,'ATG Start Percentage':atg_P,'GTG Start Percentage':gtg_P,'TTG Start Percentage':ttg_P,
                                           'ATT Start Percentage':att_P,'CTG Start Percentage':ctg_P,'Other Start Codon Percentage':other_Start_P,'TAG Stop Percentage':tag_P,'TAA Stop Percentage':taa_P,
                                           'TGA Stop Percentage':tga_P,'Other Stop Codon Percentage':other_Stop_P,'True Positive':TP,'False Positive':FP,'False Negative':FN,'Precision':precision,
                                           'Recall':recall,'False Discovery Rate':false_Discovery_Rate,'Nucleotide True Positive':NT_TP,'Nucleotide False Positive':NT_FP,'Nucleotide True Negative':NT_TN,
                                           'Nucleotide False Negative':NT_FN,'Nucleotide Precision':NT_Precision,'Nucleotide Recall':NT_Recall,'Nucleotide False Discovery Rate':NT_False_Discovery_Rate,
                                           'ORF Nucleotide Coverage of Genome':orf_Coverage_Genome,'Matched ORF Nucleotide Coverage of Genome':matched_ORF_Coverage_Genome})


    rep_Metrics = collections.OrderedDict({'Percentage of Genes Detected':genes_Detected_Percentage,'Percentage of ORFs that Detected a Gene':matched_ORF_Percentage,'Percent Difference of All ORFs':ORFs_Difference,
                                           'Median Length Difference':median_Length_Difference,'Percentage of Perfect Matches':perfect_Matches_Percentage,'Median Start Difference of Matched ORFs':median_Start_Difference,
                                           'Median Stop Difference of Matched ORFs':median_Stop_Difference,'Percent Difference of Overlapping ORFs':overlap_Difference,'Number of Short-ORFs':num_ORF_Short,
                                           'Precision':precision, 'Recall':recall,'False Discovery Rate':false_Discovery_Rate})

    #quick fix  - nan = 0

    for key,value in all_Metrics.items():
        if 'nan' == value:
            all_Metrics[key] = 'N/A'



    return all_Metrics, rep_Metrics, start_Difference, stop_Difference, other_Starts, other_Stops, comp.genes_Undetected, comp.unmatched_ORFs, Missed_Gene_Metrics, unmatched_orf_metrics,gene_Coverage_Genome, comp.multi_Matched_ORFs
