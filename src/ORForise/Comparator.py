import numpy as np

try:
    from utils import *
except ImportError:
    from ORForise.utils import *


class comparator:  # Class to hold global-type variables
    def __init__(self, perfect_Starts=0, perfect_Stops=0, genome_Seq='',
                 genome_Seq_Rev='',
                 genome_Size=0, correct_Frame_Number=0, extended_Start=0,
                 extended_Stop=0, extended_CDS=0, perfect_Matches=collections.OrderedDict(),
                 matched_ORFs=collections.OrderedDict(), multi_Matched_ORFs=collections.defaultdict(list),
                 unmatched_ORFs=collections.OrderedDict(), genes_Detected=collections.OrderedDict(),
                 genes_Undetected=collections.OrderedDict(),
                 out_Of_Frame_ORFs=collections.OrderedDict(), start_Difference=[], stop_Difference=[],
                 orf_Lengths=[], gene_Lengths=[], gene_Pos_Olap=[], gene_Neg_Olap=[], orf_Pos_Olap=[], orf_Neg_Olap=[],
                 m_ORF_Pos_Olap=[], m_ORF_Neg_Olap=[], gene_GC=[],
                 orf_GC=[], m_ORF_GC=[], gene_Short=[], orf_Short=[], m_ORF_Short=[], pos_Strand=0, neg_Strand=0,
                 partial_Hits=collections.OrderedDict()):
        self.perfect_Starts, self.perfect_Stops, self.genome_Seq, self.genome_Seq_Rev, self.genome_Size, self.correct_Frame_Number, self.extended_Start, self.extended_Stop, self.extended_CDS, \
        self.perfect_Matches, self.matched_ORFs, self.multi_Matched_ORFs, self.unmatched_ORFs, self.genes_Detected, self.genes_Undetected, self.out_Of_Frame_ORFs, self.start_Difference, \
        self.stop_Difference, self.orf_Lengths, self.gene_Lengths, self.gene_Pos_Olap, \
        self.gene_Neg_Olap, self.orf_Pos_Olap, self.orf_Neg_Olap, self.m_ORF_Pos_Olap, self.m_ORF_Neg_Olap, self.gene_GC, self.orf_GC, self.m_ORF_GC, self.gene_Short, self.orf_Short, self.m_ORF_Short, self.pos_Strand, \
        self.neg_Strand, self.partial_Hits = perfect_Starts, perfect_Stops, genome_Seq, genome_Seq_Rev, \
                                             genome_Size, correct_Frame_Number, extended_Start, extended_Stop, extended_CDS, perfect_Matches, matched_ORFs, multi_Matched_ORFs, unmatched_ORFs, genes_Detected, genes_Undetected, out_Of_Frame_ORFs, start_Difference, stop_Difference, orf_Lengths, \
                                             gene_Lengths, gene_Pos_Olap, gene_Neg_Olap, orf_Pos_Olap, orf_Neg_Olap, m_ORF_Pos_Olap, m_ORF_Neg_Olap, gene_GC, orf_GC, m_ORF_GC, gene_Short, orf_Short, m_ORF_Short, pos_Strand, neg_Strand, partial_Hits


comp = comparator()


# Not needed
# def keyshift(dictionary, key, diff):
#     if key in dictionary:
#         token = object()
#         keys = [token]*(diff*-1) + dictionary + [token]*diff
#         newkey = keys[keys.index(key)+diff]
#         if newkey is token:
#             print (None)
#         else:
#             to_return = dictionary[newkey].split(',')
#             to_return = to_return[0]+'_'+to_return[1]+'_'+to_return[2]
#             return to_return
#     else:
#         print ('Key not found')


def nuc_Count(start, stop, strand):  # Gets correct seq then returns GC
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
            t += 1
        elif "N" in i:
            n += 1
    gc_content = (g + c) * 100 / (a + t + g + c + n)
    return gc_content


def orf_Unmatched(o_Start, o_Stop, o_Strand):
    if o_Strand == '-':
        r_Start = comp.genome_Size - o_Stop
        r_Stop = comp.genome_Size - o_Start
        Unmatched_ORF = str(o_Start) + ',' + str(o_Stop) + ',' + o_Strand + ',' + comp.genome_Seq_Rev[
                                                                                  r_Start:r_Start + 3] + ',' + comp.genome_Seq_Rev[
                                                                                                               r_Stop - 2:r_Stop + 1]
        seq = (comp.genome_Seq_Rev[r_Start:r_Stop + 1])
        comp.unmatched_ORFs.update({Unmatched_ORF: seq})
    elif o_Strand == '+':
        Unmatched_ORF = str(o_Start) + ',' + str(o_Stop) + ',' + o_Strand + ',' + comp.genome_Seq[
                                                                                  o_Start - 1:o_Start + 2] + ',' + comp.genome_Seq[
                                                                                                                   o_Stop - 3:o_Stop]
        seq = (comp.genome_Seq[o_Start - 1:o_Stop])
        comp.unmatched_ORFs.update({Unmatched_ORF: seq})


def genes_Unmatched(g_Start, g_Stop, g_Strand):
    if g_Strand == '-':
        r_Start = comp.genome_Size - g_Stop
        r_Stop = comp.genome_Size - g_Start
        missed_Gene = str(g_Start) + ',' + str(g_Stop) + ',' + g_Strand + ',' + comp.genome_Seq_Rev[
                                                                                r_Start:r_Start + 3] + ',' + comp.genome_Seq_Rev[
                                                                                                             r_Stop - 2:r_Stop + 1]
        genSeq = (comp.genome_Seq_Rev[r_Start:r_Stop + 1])
        comp.genes_Undetected.update({missed_Gene: genSeq})
    elif g_Strand == '+':
        missed_Gene = str(g_Start) + ',' + str(g_Stop) + ',' + g_Strand + ',' + comp.genome_Seq[
                                                                                g_Start - 1:g_Start + 2] + ',' + comp.genome_Seq[
                                                                                                                 g_Stop - 3:g_Stop]
        genSeq = (comp.genome_Seq[g_Start - 1:g_Stop])
        comp.genes_Undetected.update({missed_Gene: genSeq})


def perfect_Matched_Genes(g_Start, g_Stop, g_Strand):
    if g_Strand == '-':
        r_Start = comp.genome_Size - g_Stop
        r_Stop = comp.genome_Size - g_Start
        perfect_Matched_Gene = str(g_Start) + ',' + str(g_Stop) + ',' + g_Strand + ',' + comp.genome_Seq_Rev[
                                                                                         r_Start:r_Start + 3] + ',' + comp.genome_Seq_Rev[
                                                                                                                      r_Stop - 2:r_Stop + 1]
        genSeq = (comp.genome_Seq_Rev[r_Start:r_Stop + 1])
        comp.perfect_Matches.update({perfect_Matched_Gene: genSeq})
    elif g_Strand == '+':
        perfect_Matched_Gene = str(g_Start) + ',' + str(g_Stop) + ',' + g_Strand + ',' + comp.genome_Seq[
                                                                                         g_Start - 1:g_Start + 2] + ',' + comp.genome_Seq[
                                                                                                                          g_Stop - 3:g_Stop]
        genSeq = (comp.genome_Seq[g_Start - 1:g_Stop])
        comp.perfect_Matches.update({perfect_Matched_Gene: genSeq})


def match_Statistics(o_Start, o_Stop, g_Start, g_Stop, g_Strand):
    comp.correct_Frame_Number += 1
    ############ Calculate prediction precision
    if '+' in g_Strand:
        comp.start_Difference.append(o_Start - g_Start)
        comp.stop_Difference.append(o_Stop - g_Stop)
        if g_Start == o_Start:
            comp.perfect_Starts += 1
        if g_Stop == o_Stop:
            comp.perfect_Stops += 1
        if o_Start < g_Start and o_Stop > g_Stop:
            comp.extended_CDS += 1
        if o_Start < g_Start:
            comp.extended_Start += 1
        if o_Stop > g_Stop:
            comp.extended_Stop += 1
    elif '-' in g_Strand:  # Negative strand genes are reversed
        comp.start_Difference.append(o_Stop - g_Stop)
        comp.stop_Difference.append(o_Start - g_Start)
        if g_Start == o_Start:
            comp.perfect_Stops += 1
        if g_Stop == o_Stop:
            comp.perfect_Starts += 1
        if o_Start < g_Start and o_Stop > g_Stop:
            comp.extended_CDS += 1
        if o_Start < g_Start:
            comp.extended_Stop += 1
        if o_Stop > g_Stop:
            comp.extended_Start += 1


def start_Codon_Count(orfs):
    atg, gtg, ttg, att, ctg, other = 0, 0, 0, 0, 0, 0
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
    atg_P = format(100 * atg / len(orfs), '.2f')
    gtg_P = format(100 * gtg / len(orfs), '.2f')
    ttg_P = format(100 * ttg / len(orfs), '.2f')
    att_P = format(100 * att / len(orfs), '.2f')
    ctg_P = format(100 * ctg / len(orfs), '.2f')
    other_Start_P = format(100 * other / len(orfs), '.2f')
    return atg_P, gtg_P, ttg_P, att_P, ctg_P, other_Start_P, other_Starts


def stop_Codon_Count(orfs):
    tag, taa, tga, other = 0, 0, 0, 0
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
    tag_p = format(100 * tag / len(orfs), '.2f')
    taa_p = format(100 * taa / len(orfs), '.2f')
    tga_p = format(100 * tga / len(orfs), '.2f')
    other_Stop_P = format(100 * other / len(orfs), '.2f')
    return tag_p, taa_p, tga_p, other_Stop_P, other_Stops


def candidate_ORF_Selection(gene_Set,
                            candidate_ORFs):  # Select ORF from candidates which is most similar to partially detected gene
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
            current_ORF_Difference = orf_Set.symmetric_difference(
                gene_Set)  # Pick least different ORF set from the Gene Set
            if len(current_ORF_Difference) > len(candidate_ORF_Difference):
                pos = c_Pos
                orf_Details = c_ORF_Details
        else:
            print("Match filtered out")
    return pos, orf_Details


def partial_Hit_Calc(g_Start, g_Stop, g_Strand, o_Start, o_Stop):
    if g_Strand == '-':
        r_G_Start = comp.genome_Size - g_Stop
        r_G_Stop = comp.genome_Size - g_Start
        r_O_Start = comp.genome_Size - o_Stop
        r_O_Stop = comp.genome_Size - o_Start
        partial = "Gene:" + str(g_Start) + '_' + str(g_Stop) + '_' + g_Strand + '_' + comp.genome_Seq_Rev[
                                                                                      r_G_Start:r_G_Start + 3] + '_' + comp.genome_Seq_Rev[
                                                                                                                       r_G_Stop - 2:r_G_Stop + 1] + ';Predicted_CDS:' + str(
            o_Start) + '_' + str(o_Stop) + '_' + g_Strand + '_' + comp.genome_Seq_Rev[
                                                                  r_O_Start:r_O_Start + 3] + '_' + comp.genome_Seq_Rev[
                                                                                                   r_O_Stop - 2:r_O_Stop + 1]
        genSeq = (comp.genome_Seq_Rev[r_G_Start:r_G_Stop + 1])
        orfSeq = (comp.genome_Seq_Rev[r_O_Start:r_O_Stop + 1])
        comp.partial_Hits.update({partial: [genSeq, orfSeq]})
    elif g_Strand == '+':
        partial = "Gene:" + str(g_Start) + '_' + str(g_Stop) + '_' + g_Strand + '_' + comp.genome_Seq[
                                                                                      g_Start - 1:g_Start + 2] + '_' + comp.genome_Seq[
                                                                                                                       g_Stop - 3:g_Stop] + ';Predicted_CDS:' + str(
            o_Start) + '_' + str(o_Stop) + '_' + g_Strand + '_' + comp.genome_Seq[
                                                                  o_Start - 1:o_Start + 2] + '_' + comp.genome_Seq[
                                                                                                   o_Stop - 3:o_Stop]
        genSeq = (comp.genome_Seq[g_Start - 1:g_Stop])
        orfSeq = (comp.genome_Seq[o_Start - 1:o_Stop])
        comp.partial_Hits.update({partial: [genSeq, orfSeq]})


def tool_comparison(ref_genes, orfs, genome, verbose):
    comp.genome_Seq = genome
    comp.genome_Seq_Rev = revCompIterative(genome)
    comp.genome_Size = len(genome)
    for gene_num, gene_details in ref_genes.items():  # Loop through each gene to compare against predicted ORFs
        g_Start = int(gene_details[0])
        g_Stop = int(gene_details[1])
        g_Strand = gene_details[2]
        g_pos = str(g_Start) + ',' + str(g_Stop)
        gene_Set = set(range(g_Start,
                             g_Stop + 1))  # Used to check Overlap of ORFs and pick best match - slow but confirms best match
        overlapping_ORFs = collections.OrderedDict()
        perfect_Match = False
        out_Frame = False
        for pos, orf_Details in orfs.items():  # Check if perfect match, if not check if match covers at least 75% of gene - Loop through ALL ORFs - SLOW
            o_Start = int(pos.split(',')[0])
            o_Stop = int(pos.split(',')[1])
            o_Strand = orf_Details[0]
            orf_Set = set(range(o_Start, o_Stop + 1))
            if o_Stop <= g_Start or o_Start >= g_Stop:  # Not caught up yet
                continue
            elif o_Start == g_Start and o_Stop == g_Stop:  # If perfect match, break and skip the rest of the ORFs
                perfect_Match = True
                break
            elif g_Start <= o_Start < g_Stop or g_Start < o_Stop < g_Stop:  # If ORF Start or Stop is between gene Start or Stop
                overlap = len(gene_Set.intersection(orf_Set))
                coverage = 100 * float(overlap) / float(len(gene_Set))
                orf_Details.append(coverage)
                if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and coverage >= MIN_COVERAGE:  # Only continue if ORF covers at least 75% of the gene and is in frame
                    overlapping_ORFs.update({pos: orf_Details})
                elif coverage >= MIN_COVERAGE:  # Not in frame / on same strand
                    comp.out_Of_Frame_ORFs.update({pos: orf_Details})
                    out_Frame = True
            elif o_Start <= g_Start and o_Stop >= g_Stop:  # If ORF extends one or both ends of the gene
                overlap = len(gene_Set.intersection(orf_Set))
                coverage = 100 * float(overlap) / float(len(gene_Set))
                orf_Details.append(coverage)
                if abs(o_Stop - g_Stop) % 3 == 0 and o_Strand == g_Strand and coverage >= MIN_COVERAGE:  # Only continue if ORF covers at least 75% of the gene and is in frame
                    overlapping_ORFs.update({pos: orf_Details})
                elif coverage >= MIN_COVERAGE:
                    comp.out_Of_Frame_ORFs.update({pos: orf_Details})
                    out_Frame = True
            else:
                if verbose == True:
                    print("Unexpected Error Finding Predicted CDSs")  # Should not happen
        # Now Check that we select the best ORF
        ### Multi_Match_ORFs Should contain All genes found by a specific ORF
        if perfect_Match == True:  # Check if the ORF is a perfect match to the Gene
            m_ORF_Details = orf_Details[:]
            m_ORF_Details.append(g_pos)
            if g_pos in comp.matched_ORFs.keys():
                previously_Covered_Gene = comp.matched_ORFs[g_pos][-1]
                comp.multi_Matched_ORFs[g_pos] += [g_pos.replace(',', '-'), previously_Covered_Gene.replace(',',
                                                                                                            '-')]  # ORF is same as gene so can use g_pos
            comp.matched_ORFs.update({g_pos: m_ORF_Details})
            comp.genes_Detected.update({str(gene_details): g_pos})
            match_Statistics(o_Start, o_Stop, g_Start, g_Stop, g_Strand)
            perfect_Matched_Genes(g_Start, g_Stop, g_Strand)
            if verbose == True:
                print('Perfect Match')
        elif perfect_Match == False and len(
                overlapping_ORFs) == 1:  # If we do not have a perfect match but 1 ORF which has passed the filtering
            orf_Pos = list(overlapping_ORFs.keys())[0]
            o_Start = int(orf_Pos.split(',')[0])
            o_Stop = int(orf_Pos.split(',')[1])
            orf_Details = overlapping_ORFs[orf_Pos]
            m_ORF_Details = orf_Details[:]
            m_ORF_Details.append(g_pos)
            if orf_Pos in comp.matched_ORFs.keys():
                try:
                    previously_Covered_Gene = comp.matched_ORFs[g_pos][-1]
                except KeyError:
                    last_key = [*comp.matched_ORFs.keys()][-1]
                    previously_Covered_Gene = comp.matched_ORFs[last_key][-1]
                comp.multi_Matched_ORFs[orf_Pos] += [g_pos.replace(',', '-'), previously_Covered_Gene.replace(',',
                                                                                                   '-')]  # ORF collects multiple gene pos'
            comp.matched_ORFs.update({orf_Pos: m_ORF_Details})
            comp.genes_Detected.update({str(gene_details): orf_Pos})
            match_Statistics(o_Start, o_Stop, g_Start, g_Stop, g_Strand)
            if verbose == True:
                print('Partial Match')
            partial_Hit_Calc(g_Start, g_Stop, g_Strand, o_Start, o_Stop)
        elif perfect_Match == False and len(
                overlapping_ORFs) >= 1:  # If we have more than 1 potential ORF match, we check to see which is the 'best' hit
            orf_Pos, orf_Details = candidate_ORF_Selection(gene_Set, overlapping_ORFs)  # Return best match
            o_Start = int(orf_Pos.split(',')[0])
            o_Stop = int(orf_Pos.split(',')[1])
            m_ORF_Details = orf_Details[:]
            m_ORF_Details.append(g_pos)
            if orf_Pos in comp.matched_ORFs.keys():
                try:
                    previously_Covered_Gene = comp.matched_ORFs[g_pos][-1]
                except KeyError:
                    last_key = [*comp.matched_ORFs.keys()][-1]
                    previously_Covered_Gene = comp.matched_ORFs[last_key][-1]
                comp.multi_Matched_ORFs[orf_Pos] += [g_pos.replace(',', '-'), previously_Covered_Gene.replace(',',
                                                                                                              '-')]  # ORF collects multiple gene pos'
            comp.matched_ORFs.update({orf_Pos: m_ORF_Details})
            comp.genes_Detected.update({str(gene_details): orf_Pos})
            match_Statistics(o_Start, o_Stop, g_Start, g_Stop, g_Strand)
            if verbose == True:
                print('There was more than 1 potential Match - Best Chosen')
            partial_Hit_Calc(g_Start, g_Stop, g_Strand, o_Start, o_Stop)
        elif out_Frame:  # Keep record of ORFs which overlap a gene but in the wrong frame
            if verbose == True:
                print("Out of Frame Predicted CDS")
            genes_Unmatched(g_Start, g_Stop, g_Strand)  #
        else:
            genes_Unmatched(g_Start, g_Stop, g_Strand)  # No hit
            if verbose == True:
                print("No Hit")
    for orf_Key in comp.matched_ORFs:  # Remove ORFs from out of frame if ORF was correctly matched to another Gene
        if orf_Key in comp.out_Of_Frame_ORFs:
            del comp.out_Of_Frame_ORFs[orf_Key]
    ######################################## ORF Lengths and Precision
    start_Difference = [x for x in comp.start_Difference if x != 0]  # Remove 0s (Perfect hits)
    stop_Difference = [x for x in comp.stop_Difference if x != 0]
    if len(start_Difference) >= 1:
        median_Start_Difference = np.median(start_Difference)
    else:
        median_Start_Difference = 'N/A'
    if len(stop_Difference) >= 1:
        median_Stop_Difference = np.median(stop_Difference)
    else:
        median_Stop_Difference = 'N/A'

    # Get Start and Stop Codon Usage
    atg_P, gtg_P, ttg_P, att_P, ctg_P, other_Start_P, other_Starts = start_Codon_Count(orfs)
    tag_P, taa_P, tga_P, other_Stop_P, other_Stops = stop_Codon_Count(orfs)
    # Count nucleotides found from ALL ORFs
    gene_Nuc_Array = np.zeros((comp.genome_Size), dtype=np.int)
    orf_Nuc_Array = np.zeros((comp.genome_Size), dtype=np.int)
    matched_ORF_Nuc_Array = np.zeros((comp.genome_Size), dtype=np.int)

    prev_Gene_Stop = 0
    prev_Gene_Overlapped = False
    for gene_Num, gene_details in ref_genes.items():  # Loop through each gene to compare against predicted ORFs
        g_Start = int(gene_details[0])
        g_Stop = int(gene_details[1])
        g_Strand = gene_details[2]
        gene_Length = (g_Stop - g_Start)
        comp.gene_Lengths.append(gene_Length)
        gene_Nuc_Array[g_Start - 1:g_Stop] = [1]  # Changing all between the two positions to 1's
        comp.gene_GC.append(nuc_Count(g_Start, g_Stop, g_Strand))
        if gene_Length <= SHORT_ORF_LENGTH:  # .utils
            comp.gene_Short.append(gene_Length)
        ### Calculate overlapping Genes -
        if prev_Gene_Stop > g_Start:
            if '+' in g_Strand:
                comp.gene_Pos_Olap.append(prev_Gene_Stop - g_Start)
            elif '-' in g_Strand:
                comp.gene_Neg_Olap.append(prev_Gene_Stop - g_Start)
            prev_Gene_Overlapped = True
        elif prev_Gene_Stop < g_Start:
            if prev_Gene_Overlapped == True:
                if '+' in g_Strand:
                    comp.gene_Pos_Olap.append(0)
                elif '-' in g_Strand:
                    comp.gene_Neg_Olap.append(0)
            prev_Gene_Overlapped = False
        prev_Gene_Stop = g_Stop
    if prev_Gene_Overlapped == True:  # If last has a prev overlap, count it
        if '+' in g_Strand:
            comp.gene_Pos_Olap.append(0)
        elif '-' in g_Strand:
            comp.gene_Neg_Olap.append(0)
    ####
    min_Gene_Length = min(comp.gene_Lengths)
    max_Gene_Length = max(comp.gene_Lengths)
    median_Gene_Length = np.median(comp.gene_Lengths)
    prev_ORF_Stop = 0
    prev_ORF_Overlapped = False
    for o_Positions, orf_Details in orfs.items():
        o_Start = int(o_Positions.split(',')[0])
        o_Stop = int(o_Positions.split(',')[1])
        o_Strand = orf_Details[0]
        # Stats just for Unmatched ORFs
        if o_Positions not in list(comp.matched_ORFs.keys()):
            orf_Unmatched(o_Start, o_Stop, o_Strand)
        # Get ORF Strand metrics:
        if o_Strand == "+":  # Get number of Positive and Negative strand ORFs
            comp.pos_Strand += 1
        elif o_Strand == "-":
            comp.neg_Strand += 1
        orf_Length = (o_Stop - o_Start)
        comp.orf_Lengths.append(orf_Length)
        orf_Nuc_Array[o_Start - 1:o_Stop] = [1]  # Changing all between the two positions to 1's
        comp.orf_GC.append(nuc_Count(o_Start, o_Stop, o_Strand))
        if orf_Length <= SHORT_ORF_LENGTH:  # .utils
            comp.orf_Short.append(orf_Length)
        ### Calculate overlapping ORFs -
        if prev_ORF_Stop > o_Start:
            if '+' in o_Strand:
                comp.orf_Pos_Olap.append(prev_ORF_Stop - o_Start)
            elif '-' in o_Strand:
                comp.orf_Neg_Olap.append(prev_ORF_Stop - o_Start)
            prev_ORF_Overlapped = True
        elif prev_ORF_Stop < o_Start:
            if prev_ORF_Overlapped == True:
                if '+' in o_Strand:
                    comp.orf_Pos_Olap.append(0)
                elif '-' in o_Strand:
                    comp.orf_Neg_Olap.append(0)
            prev_ORF_Overlapped = False
        prev_ORF_Stop = o_Stop
        if prev_ORF_Overlapped == True:  # If last has a prev overlap, count it
            if '+' in o_Strand:
                comp.orf_Pos_Olap.append(0)
            elif '-' in o_Strand:
                comp.orf_Neg_Olap.append(0)

    # Nucleotide Coverage calculated from ORFs matching a gene only
    matched_Prev_ORF_Stop = 0
    matched_Prev_ORF_Overlapped = False
    for mo_Positions, m_ORF_Details in comp.matched_ORFs.items():
        mo_Start = int(mo_Positions.split(',')[0])
        mo_Stop = int(mo_Positions.split(',')[1])
        mo_Strand = m_ORF_Details[0]
        mo_Length = (mo_Stop - mo_Start)
        matched_ORF_Nuc_Array[mo_Start - 1:mo_Stop] = [1]  # This is the complete matched orf not the matched orf bits

        comp.m_ORF_GC.append(nuc_Count(mo_Start, mo_Stop, mo_Strand))
        if mo_Length <= SHORT_ORF_LENGTH:  # .utils
            comp.m_ORF_Short.append(mo_Length)
        ### Calculate overlapping Matched ORFs -
        if matched_Prev_ORF_Stop > mo_Start:
            if '+' in mo_Strand:
                comp.m_ORF_Pos_Olap.append(matched_Prev_ORF_Stop - mo_Start)
            elif '-' in mo_Strand:
                comp.m_ORF_Neg_Olap.append(matched_Prev_ORF_Stop - mo_Start)
            matched_Prev_ORF_Overlapped = True
        elif matched_Prev_ORF_Stop < mo_Start:
            if matched_Prev_ORF_Overlapped == True:
                if '+' in mo_Strand:
                    comp.m_ORF_Pos_Olap.append(0)
                elif '-' in mo_Strand:
                    comp.m_ORF_Neg_Olap.append(0)
            matched_Prev_ORF_Overlapped = False
        matched_Prev_ORF_Stop = mo_Stop
    if matched_Prev_ORF_Overlapped == True:  # If last has a prev overlap, count it
        if '+' in mo_Strand:
            comp.m_ORF_Pos_Olap.append(0)
        elif '-' in mo_Strand:
            comp.m_ORF_Neg_Olap.append(0)
    ####
    gene_Coverage_Genome = format(100 * np.count_nonzero(gene_Nuc_Array) / comp.genome_Size, '.2f')
    orf_Coverage_Genome = format(100 * np.count_nonzero(orf_Nuc_Array) / comp.genome_Size, '.2f')
    matched_ORF_Coverage_Genome = format(100 * np.count_nonzero(matched_ORF_Nuc_Array) / comp.genome_Size,
                                         '.2f')  # This gets the nts which are in matched ORFs - Check below
    # matched_ORF_Nuc_AND_Gene = np.logical_and(matched_ORF_Nuc_Array,gene_Nuc_Array) + [0 for i in range(len(gene_Nuc_Array))] # This gets the nts which are in both matched ORFs and detected genes
    # matched_ORF_Coverage_Genome = format(100 * np.count_nonzero(matched_ORF_Nuc_AND_Gene) / comp.genome_Size,'.2f')

    # gene and orf nucleotide Intersection
    gene_ORF_Nuc_Intersection = np.count_nonzero(gene_Nuc_Array & orf_Nuc_Array)
    # not gene but orf nucleotides
    not_Gene_Nuc_Array = np.logical_not(gene_Nuc_Array) + [0 for i in range(
        len(gene_Nuc_Array))]  # End part to keep array as 1,0 not T,F
    not_Gene_Nuc_And_ORF_Count = np.count_nonzero(not_Gene_Nuc_Array & orf_Nuc_Array)
    # not orf nucleotides but gene
    not_ORF_Nuc_Array = np.logical_not(orf_Nuc_Array) + [0 for i in range(
        len(orf_Nuc_Array))]  # End part to keep array as 1,0 not T,F
    not_ORF_Nuc_And_Gene_Count = np.count_nonzero(not_ORF_Nuc_Array & gene_Nuc_Array)
    # not gene or orf nucleotides
    not_Gene_Nuc_Not_ORF_Nuc_Count = np.count_nonzero(not_Gene_Nuc_Array & not_ORF_Nuc_Array)
    # Nucleotide 'accuracy' - Normalised by number of nucelotides annotated by a gene
    NT_TP = format(gene_ORF_Nuc_Intersection / np.count_nonzero(gene_Nuc_Array), '.2f')
    NT_FP = format(not_Gene_Nuc_And_ORF_Count / np.count_nonzero(not_Gene_Nuc_Array), '.2f')
    NT_FN = format(not_ORF_Nuc_And_Gene_Count / np.count_nonzero(gene_Nuc_Array), '.2f')
    NT_TN = format(not_Gene_Nuc_Not_ORF_Nuc_Count / np.count_nonzero(not_Gene_Nuc_Array), '.2f')
    NT_Precision = format(gene_ORF_Nuc_Intersection / (gene_ORF_Nuc_Intersection + not_Gene_Nuc_And_ORF_Count), '.2f')
    NT_Recall = format(gene_ORF_Nuc_Intersection / (gene_ORF_Nuc_Intersection + not_ORF_Nuc_And_Gene_Count), '.2f')
    NT_False_Discovery_Rate = format(
        not_Gene_Nuc_And_ORF_Count / (not_Gene_Nuc_And_ORF_Count + gene_ORF_Nuc_Intersection), '.2f')
    ################################# Precision and Recall of whole ORFs and Genes
    TP = format(len(comp.genes_Detected) / len(ref_genes), '.2f')
    FP = format(len(comp.unmatched_ORFs) / len(ref_genes), '.2f')
    FN = format(len(comp.genes_Undetected) / len(ref_genes), '.2f')
    #################################################### Need a better way to handle 'no hits/ORFs'
    try:
        precision = format(float(TP) / (float(TP) + float(FP)), '.2f')
    except ZeroDivisionError:
        precision = format(0.00, '.2f')
    try:
        recall = format(float(TP) / (float(TP) + float(FN)), '.2f')
    except ZeroDivisionError:
        recall = format(0.00, '.2f')
    try:
        false_Discovery_Rate = format(float(FP) / (float(FP) + float(TP)), '.2f')
    except ZeroDivisionError:
        false_Discovery_Rate = 'N/A'
    min_ORF_Length = min(comp.orf_Lengths)
    max_ORF_Length = max(comp.orf_Lengths)
    median_ORF_Length = np.median(comp.orf_Lengths)

    ##########################################################################
    # Metrics - There are numerous cases where certain metric calculations may return a ZeroDivError.
    ORFs_Difference = format(100 * (len(orfs) - len(ref_genes)) / len(ref_genes), '.2f')  # Difference off +/-
    genes_Detected_Percentage = format(100 * (len(comp.genes_Detected) / len(ref_genes)), '.2f')
    matched_ORF_Percentage = format(100 * (len(comp.matched_ORFs) / len(orfs)), '.2f')
    all_ORF_Olap = (comp.orf_Pos_Olap + comp.orf_Neg_Olap)  # Combine pos and neg strand overlaps
    matched_ORF_Olap = (comp.m_ORF_Pos_Olap + comp.m_ORF_Neg_Olap)
    all_Gene_Olap = (comp.gene_Pos_Olap + comp.gene_Neg_Olap)

    if all_ORF_Olap:  # If no overlapping ORFs
        try:
            overlap_Difference = format(100 * (len(all_ORF_Olap) - len(all_Gene_Olap)) / len(all_Gene_Olap), '.2f')
            matched_Overlap_Difference = format(100 * (len(matched_ORF_Olap) - len(all_Gene_Olap)) / len(all_Gene_Olap),
                                                '.2f')
        except ZeroDivisionError:
            overlap_Difference = 'N/A'
            matched_Overlap_Difference = 'N/A'
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
        overlap_Difference = 'N/A'
        matched_Overlap_Difference = 'N/A'
        num_All_ORF_Olap = 0
        max_Matched_ORF_Olap = 'N/A'
        max_All_ORF_Olap = 'N/A'
        median_ORF_Overlap = 'N/A'
        matched_Median_ORF_Overlap = 'N/A'
    if len(matched_ORF_Olap) == 0:  # -100.00 is not informative
        matched_Overlap_Difference = 'N/A'

    # Need to NA everything

    if comp.orf_Short and comp.gene_Short:  # IF Short-ORFs/Genes
        short_ORF_Difference = format(100 * (len(comp.orf_Short) - len(comp.gene_Short)) / len(comp.gene_Short), '.2f')
        matched_Short_ORF_Difference = format(
            100 * (len(comp.m_ORF_Short) - len(comp.gene_Short)) / len(comp.gene_Short), '.2f')
        num_ORF_Short = len(comp.orf_Short)
        num_Matched_ORF_Short = len(comp.m_ORF_Short)
    elif comp.orf_Short:  # If only Short-ORFs
        num_ORF_Short = len(comp.orf_Short)
        num_Matched_ORF_Short = 'N/A'
        short_ORF_Difference = (num_ORF_Short * 100)
        matched_Short_ORF_Difference = 'N/A'
    else:  # If only Short-Genes and Undetected StORFs
        comp.gene_Short
        short_ORF_Difference = 'N/A'
        matched_Short_ORF_Difference = 'N/A'
        num_ORF_Short = 0
        num_Matched_ORF_Short = 'N/A'
    if num_Matched_ORF_Short == 0:  # -100.00 is not informative
        matched_Short_ORF_Difference = 'N/A'

    median_Length_Difference = format(100 * (median_ORF_Length - median_Gene_Length) / median_Gene_Length, '.2f')
    min_Length_Difference = format(100 * (min_ORF_Length - min_Gene_Length) / min_Gene_Length, '.2f')
    max_Length_Difference = format(100 * (max_ORF_Length - max_Gene_Length) / max_Gene_Length, '.2f')
    pos_Strand_Percentage = format(comp.pos_Strand / len(orfs), '.2f')
    neg_Strand_Percentage = format(comp.neg_Strand / len(orfs), '.2f')
    median_ORF_GC = np.median(comp.orf_GC)
    matched_Median_ORF_GC = np.median(comp.m_ORF_GC)
    median_Gene_GC = np.median(comp.gene_GC)
    median_GC_Difference = format(100 * (float(median_ORF_GC) - float(median_Gene_GC)) / float(median_Gene_GC), '.2f')
    matched_Median_GC_Difference = format(
        100 * (float(matched_Median_ORF_GC) - float(median_Gene_GC)) / float(median_Gene_GC), '.2f')

    if comp.matched_ORFs:  # No ORFs detected a gene
        extended_CDS_Percentage = format(100 * comp.extended_CDS / len(comp.matched_ORFs), '.2f')
        extended_Start_Percentage = format(100 * comp.extended_Start / len(comp.matched_ORFs), '.2f')
        extended_Stop_Percentage = format(100 * comp.extended_Stop / len(comp.matched_ORFs), '.2f')
        perfect_Matches_Percentage = format(100 * len(comp.perfect_Matches) / len(comp.matched_ORFs), '.2f')
        perfect_Starts_Percentage = format(100 * comp.perfect_Starts / len(comp.matched_ORFs), '.2f')
        perfect_Stops_Percentage = format(100 * comp.perfect_Stops / len(comp.matched_ORFs), '.2f')
    else:
        # correct_Frame_Percentage = 0
        extended_CDS_Percentage = format(0.00, '.2f')
        extended_Start_Percentage = format(0.00, '.2f')
        extended_Stop_Percentage = format(0.00, '.2f')
        perfect_Matches_Percentage = format(0.00, '.2f')
        perfect_Starts_Percentage = format(0.00, '.2f')
        perfect_Stops_Percentage = format(0.00, '.2f')
    ################### Missed Genes  Metrics:
    if comp.genes_Undetected:
        mg_Starts = []
        mg_Stops = []
        mg_Lengths = []
        mg_Strands = []
        for mg, seq in comp.genes_Undetected.items():
            mg = mg.split(',')
            mg_Starts.append(mg[3])
            mg_Stops.append(mg[4])
            mg_Strands.append(mg[2])
            mg_Lengths.append(int(mg[1]) - int(mg[0]))

        mg_ATG = 100 * mg_Starts.count('ATG') / len(comp.genes_Undetected)
        mg_GTG = 100 * mg_Starts.count('GTG') / len(comp.genes_Undetected)
        mg_TTG = 100 * mg_Starts.count('TTG') / len(comp.genes_Undetected)
        mg_ATT = 100 * mg_Starts.count('ATT') / len(comp.genes_Undetected)
        mg_CTG = 100 * mg_Starts.count('CTG') / len(comp.genes_Undetected)
        mg_O_Start = 100 - (mg_ATG + mg_GTG + mg_TTG + mg_ATT + mg_CTG)
        mg_TGA = 100 * mg_Stops.count('TGA') / len(comp.genes_Undetected)
        mg_TAA = 100 * mg_Stops.count('TAA') / len(comp.genes_Undetected)
        mg_TAG = 100 * mg_Stops.count('TAG') / len(comp.genes_Undetected)
        mg_O_Stop = 100 - (mg_TGA + mg_TAA + mg_TAG)
        median_mg_Len = np.median(mg_Lengths)
        mg_Pos = mg_Strands.count('+')
        mg_Neg = mg_Strands.count('-')
        undetected_Gene_Metrics = (
        format(mg_ATG, '.2f'), format(mg_GTG, '.2f'), format(mg_TTG, '.2f'), format(mg_ATT, '.2f'),
        format(mg_CTG, '.2f'), format(mg_O_Start, '.2f'), format(mg_TGA, '.2f'), format(mg_TAA, '.2f'),
        format(mg_TAG, '.2f'), format(mg_O_Stop, '.2f'), format(median_mg_Len, '.2f'), mg_Pos, mg_Neg)
    else:
        undetected_Gene_Metrics = ''
    #################### Unmathced ORF Metrics:
    if comp.unmatched_ORFs:
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
        uo_O_Stop = 100 - (uo_TGA + uo_TAA + uo_TAG)
        # uo_O_Stop = 100 * uo_O_Stop  / len(comp.unmatched_ORFs) ########WHY?
        median_uo_Len = np.median(uo_Lengths)
        uo_Pos = uo_Strands.count('+')
        uo_Neg = uo_Strands.count('-')
        unmatched_ORF_Metrics = (
        format(uo_ATG, '.2f'), format(uo_GTG, '.2f'), format(uo_TTG, '.2f'), format(uo_ATT, '.2f'),
        format(uo_CTG, '.2f'), format(uo_O_Start, '.2f'), format(uo_TGA, '.2f'), format(uo_TAA, '.2f'),
        format(uo_TAG, '.2f'), format(uo_O_Stop, '.2f'), format(median_uo_Len, '.2f'), uo_Pos, uo_Neg)
    else:
        unmatched_ORF_Metrics = ''
    #################################

    all_Metrics = collections.OrderedDict(
        {'Number_of_ORFs': len(orfs), 'Percent_Difference_of_All_ORFs': ORFs_Difference,
         'Number_of_ORFs_that_Detected_a_Gene': len(comp.matched_ORFs),
         'Percentage_of_ORFs_that_Detected_a_Gene': matched_ORF_Percentage,
         'Number_of_Genes_Detected': len(comp.genes_Detected),
         'Percentage_of_Genes_Detected': genes_Detected_Percentage, 'Median_Length_of_All_ORFs': median_ORF_Length,
         'Median_Length_Difference': median_Length_Difference,
         'Minimum_Length_of_All_ORFs': min_ORF_Length, 'Minimum_Length_Difference': min_Length_Difference,
         'Maximum_Length_of_All_ORFs': max_ORF_Length, 'Maximum_Length_Difference': max_Length_Difference,
         'Median_GC_content_of_All_ORFs': format(median_ORF_GC, '.2f'),
         'Percent_Difference_of_All_ORFs_Median_GC': median_GC_Difference,
         'Median_GC_content_of_Matched_ORFs': format(matched_Median_ORF_GC, '.2f'),
         'Percent_Difference_of_Matched_ORF_GC': matched_Median_GC_Difference,
         'Number_of_ORFs_which_Overlap_Another_ORF': num_All_ORF_Olap,
         'Percent_Difference_of_Overlapping_ORFs': overlap_Difference,
         'Maximum_ORF_Overlap': max_All_ORF_Olap, 'Median_ORF_Overlap': median_ORF_Overlap,
         'Number_of_Matched_ORFs_Overlapping_Another_ORF': len(matched_ORF_Olap),
         'Percentage_Difference_of_Matched_Overlapping_CDSs': matched_Overlap_Difference,
         'Maximum_Matched_ORF_Overlap': max_Matched_ORF_Olap, 'Median_Matched_ORF_Overlap': matched_Median_ORF_Overlap,
         'Number_of_Short-ORFs': num_ORF_Short, 'Percent_Difference_of_Short-ORFs': short_ORF_Difference,
         'Number_of_Short-Matched-ORFs': num_Matched_ORF_Short,
         'Percent_Difference_of_Short-Matched-ORFs': matched_Short_ORF_Difference,
         'Number_of_Perfect_Matches': len(comp.perfect_Matches),
         'Percentage_of_Perfect_Matches': perfect_Matches_Percentage,
         'Number_of_Perfect_Starts': comp.perfect_Starts, 'Percentage_of_Perfect_Starts': perfect_Starts_Percentage,
         'Number_of_Perfect_Stops': comp.perfect_Stops, 'Percentage_of_Perfect_Stops': perfect_Stops_Percentage,
         'Number_of_Out_of_Frame_ORFs': len(comp.out_Of_Frame_ORFs),
         'Number_of_Matched_ORFs_Extending_a_Coding_Region': comp.extended_CDS,
         'Percentage_of_Matched_ORFs_Extending_a_Coding_Region': extended_CDS_Percentage,
         'Number_of_Matched_ORFs_Extending_Start_Region': comp.extended_Start,
         'Percentage_of_Matched_ORFs_Extending_Start_Region': extended_Start_Percentage,
         'Number_of_Matched_ORFs_Extending_Stop_Region': comp.extended_Stop,
         'Percentage_of_Matched_ORFs_Extending_Stop_Region': extended_Stop_Percentage,
         'Number_of_All_ORFs_on_Positive_Strand': comp.pos_Strand,
         'Percentage_of_All_ORFs_on_Positive_Strand': pos_Strand_Percentage,
         'Number_of_All_ORFs_on_Negative_Strand': comp.neg_Strand,
         'Percentage_of_All_ORFs_on_Negative_Strand': neg_Strand_Percentage,
         'Median_Start_Difference_of_Matched_ORFs': median_Start_Difference,
         'Median_Stop_Difference_of_Matched_ORFs': median_Stop_Difference, 'ATG_Start_Percentage': atg_P,
         'GTG_Start_Percentage': gtg_P, 'TTG_Start_Percentage': ttg_P,
         'ATT_Start_Percentage': att_P, 'CTG_Start_Percentage': ctg_P, 'Other_Start_Codon_Percentage': other_Start_P,
         'TAG_Stop_Percentage': tag_P, 'TAA_Stop_Percentage': taa_P,
         'TGA_Stop_Percentage': tga_P, 'Other_Stop_Codon_Percentage': other_Stop_P, 'True_Positive': TP,
         'False_Positive': FP, 'False_Negative': FN, 'Precision': precision,
         'Recall': recall, 'False_Discovery_Rate': false_Discovery_Rate, 'Nucleotide_True_Positive': NT_TP,
         'Nucleotide_False_Positive': NT_FP, 'Nucleotide_True_Negative': NT_TN,
         'Nucleotide_False_Negative': NT_FN, 'Nucleotide_Precision': NT_Precision, 'Nucleotide_Recall': NT_Recall,
         'Nucleotide_False_Discovery_Rate': NT_False_Discovery_Rate,
         'ORF_Nucleotide_Coverage_of_Genome': orf_Coverage_Genome,
         'Matched_ORF_Nucleotide_Coverage_of_Genome': matched_ORF_Coverage_Genome})

    rep_Metrics = collections.OrderedDict(
        {'Percentage_of_Genes_Detected': genes_Detected_Percentage,
         'Percentage_of_ORFs_that_Detected_a_Gene': matched_ORF_Percentage,
         'Percent_Difference_of_All_ORFs': ORFs_Difference,
         'Median_Length_Difference': median_Length_Difference,
         'Percentage_of_Perfect_Matches': perfect_Matches_Percentage,
         'Median_Start_Difference_of_Matched_ORFs': median_Start_Difference,
         'Median_Stop_Difference_of_Matched_ORFs': median_Stop_Difference,
         'Percentage_Difference_of_Matched_Overlapping_CDSs': matched_Overlap_Difference,
         'Percent_Difference_of_Short-Matched-ORFs': matched_Short_ORF_Difference,
         'Precision': precision,
         'Recall': recall,
         'False_Discovery_Rate': false_Discovery_Rate})

    # To account for unbalanced data
    for m_key, m_value in all_Metrics.items():
        if 'nan' == m_value:
            all_Metrics[m_key] = 'N/A'

    return all_Metrics, rep_Metrics, start_Difference, stop_Difference, other_Starts, other_Stops, comp.perfect_Matches, comp.genes_Undetected, comp.unmatched_ORFs, undetected_Gene_Metrics, unmatched_ORF_Metrics, orf_Coverage_Genome, matched_ORF_Coverage_Genome, gene_Coverage_Genome, comp.multi_Matched_ORFs, comp.partial_Hits
