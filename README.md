# ORForise - Prokaryote Genome Annotation Analysis and Comparison Platform
## Published in Bioinformatics :   https://academic.oup.com/bioinformatics/article/38/5/1198/6454948
### Platform for analysing and comparing Prokaryote CoDing Sequence (CDS) Gene Predictions. 
### Novel genome annotations can be compared to a provided reference annotation from Ensembl and predictions from other tools (or any given GFF annotation) .

# Requirements and Installation:

### The ORForise platform is written in Python (3.6-3.9) and only requires the NumPy library (should be installed automatically by pip when installing ORForise) which is standard in most base installations of Python3.

## Intallation:

### The ORForise platform is available via the pip Python package manager ```pip3 install ORForise```. 
### Consider using '--no-cache-dir' with pip to ensure the download of the newest version of the package.

## Required Files:

To run, you need:
* Input Genome FASTA and corresponding GFF file (or CDS predictions with the annotated genes for the genome you want to use as reference in one of the tool output formats listed below).
* A prediction output from one of the compatible tools for the same genome.

### How to add your own Genome:

Corresponding FASTA and GFF files must be provided for the genome the analysis is to be performed on, including the corresponding output of any tools to compare.

### How to add your own tool:

If the new tool reports its predictions in GFF you can present ORForise with "GFF" for either the reference ```-rt``` or prediction ```-t``` option.
If the tool uses another non-standard format, a request can be made to add it as an option via GitHub.


### Testing:
Precomputed testing and data which includes example input and output files for all tools presented below is available in the `~ORForise/Testing` directory of the GitHub repository. 
Example output files from ```Annotation-Compare```, ```GFF-Adder``` and ```GFF-Intersector``` are made available to validate installation.


## CDS Prediction Analysis:

### Use-cases: (Running if via pip)

For Help: ```Annotation-Compare -h ```

```python
Thank you for using ORForise
Please report any issues to: https://github.com/NickJD/ORForise/issues
#####
usage: Annotation_Compare.py [-h] -dna GENOME_DNA -ref REFERENCE_ANNOTATION -t TOOL -tp TOOL_PREDICTION
                             [-rt REFERENCE_TOOL] [-o OUTNAME] [-v {True,False}]

ORForise v1.4.1: Annotatione-Compare Run Parameters.

Required Arguments:
  -dna GENOME_DNA       Genome DNA file (.fa) which both annotations are based on
  -ref REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -t TOOL               Which tool to analyse? (Prodigal)
  -tp TOOL_PREDICTION   Tool genome prediction file (.gff) - Different Tool Parameters are compared individually via
                        separate files

Optional Arguments:
  -rt REFERENCE_TOOL    What type of Annotation to compare to? -- Leave blank for Ensembl reference- Provide tool
                        name to compare output from two tools

Output:
  -o OUTNAME            Define full output filename (format is CSV) - If not provided, summary will be printed to
                        std-out

Misc:
  -v {True,False}       Default - False: Print out runtime status
```

## Compare a novel genome annotation to an Ensembl annotation:

Genome annotation is a difficult process, even for Prokaryotes. ORForise allows the direct and systematic analysis of
a novel CDS prediction from a wide selection of tools to a reference Genome Annotation, such as those provided by
Ensembl Bacteria.

#### Example: Installation through pip will allow user to call the programs directly from the ORForise package.
```python
Annotation-Compare -dna ~/Testing/Myco.fa -ref ~/Testing/Myco.gff -t Prodigal -tp ~/Testing/Prodigal_Myco.gff
```
### Compare different novel annotations with each other on a single Genome:

If a reference Genome Annotation is not available or a direct comparison between two or more tools is wanted,
ORForise can be used as the example below.

## Aggregate CDS Prediction Analysis:

### Use-cases: (Running if via pip)

For Help: ```Aggregate-Compare -h ```

```python
Thank you for using ORForise
Please report any issues to: https://github.com/NickJD/ORForise/issues
#####
usage: Aggregate_Compare.py [-h] -dna GENOME_DNA -t TOOLS -tp TOOL_PREDICTIONS -ref REFERENCE_ANNOTATION
                            [-rt REFERENCE_TOOL] [-o OUTNAME] [-v {True,False}]

ORForise v1.4.1: Aggregate-Compare Run Parameters.

Required Arguments:
  -dna GENOME_DNA       Genome DNA file (.fa) which both annotations are based on
  -t TOOLS              Which tools to analyse? (Prodigal,GeneMarkS)
  -tp TOOL_PREDICTIONS  Tool genome prediction file (.gff) - Providefile locations for each tool comma separated
  -ref REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?

Optional Arguments:
  -rt REFERENCE_TOOL    What type of Annotation to compare to? -- Leave blank for Ensembl reference- Provide tool
                        name to compare output from two tools

Output:
  -o OUTNAME            Define full output filename (format is CSV) - If not provided, summary will be printed to
                        std-out

Misc:
  -v {True,False}       Default - False: Print out runtime status
```

#### Example: 
```python
Aggregate-Compare -ref ~/Testing/Myco.gff -dna ~/Testing/Myco.fa -t Prodigal,TransDecoder,GeneMark_S_2 -tp ~/Testing/Prodigal_Myco.gff,~/Testing/TransDecoder_Myco.gff,~/Testing/GeneMark_S_2_Myco.gff
```
This will compare the Aggregate the predictions of Prodigal, TransDecoder and GLIMMER 3 against the Mycoplasma reference annotation provided by
Ensembl Bacteria.

## Annotation Comparison Output - The output format is the same for Annotation_Compare and Aggregate_Compare:
### Print to screen example - Prodigal prediction compared to Ensembl Bacteria reference annotation of *Escherichia coli*:
```bash
Annotation-Compare.py  -ref ./Testing/Myco.gff -dna ./Testing/Myco.fa -t Prodigal -tp ./Testing/Prodigal_Myco.gff
Genome Used: Myco
Reference Used: Testing/Myco.gff
Tool Compared: Prodigal
Perfect Matches:128[476] -26.89%
Partial Matches:62[476] - 13.03%
Missed Genes:286[476] - 60.08%
Complete
```

``` bash
Aggregate-Compare -ref ./Testing/Myco.gff -dna ./Testing/Myco.fa -t Prodigal,TransDecoder,GeneMark_S_2 -tp ./Testing/Prodigal_Myco.gff,./Testing/TransDecoder_Myco.gff,./Testing/GeneMark_S_2_Myco.gff
Prodigal
TransDecoder
GeneMark_S_2
Match filtered out
Match filtered out
Match filtered out
Match filtered out
Match filtered out
Match filtered out
Genome Used: Myco
Reference Used: ./Testing/Myco.gff
Tools Compared: Prodigal,TransDecoder,GeneMark_S_2
Perfect Matches:132[476]
Partial Matches:58[476]
Missed Genes:286[476]
```

This is the default output of the comparison tools. 

### '-o' Example output to CSV file - Prodigal prediction compared to Ensembl Bacteria reference annotation of *Escherichia coli*:
The output is designed to be human-readable and interpretable by the included 'ORForise_Analysis' scripts. 
The example below presents the 12 'Representative' and 72 'All' Metrics but only shows one entry for each of the induvidual prediction reports (Perfect_Match_Genes,Partial_Match_Genes,Missed_Genes,Predicted_CDS_Without_Corresponding_Gene_in_Reference,Predicted_CDSs_Which_Detected_more_than_one_Gene).

```csv
Representative_Metrics:
Percentage_of_Genes_Detected,Percentage_of_ORFs_that_Detected_a_Gene,Percent_Difference_of_All_ORFs,Median_Length_Difference,Percentage_of_Perfect_Matches,Median_Start_Difference_of_Matched_ORFs,Median_Stop_Difference_of_Matched_ORFs,Percentage_Difference_of_Matched_Overlapping_CDSs,Percent_Difference_of_Short-Matched-ORFs,Precision,Recall,False_Discovery_Rate
39.92,19.10,109.03,-62.17,67.37,67.5,-85.5,-83.71,-17.39,0.19,0.40,0.81
All_Metrics:
Number_of_ORFs,Percent_Difference_of_All_ORFs,Number_of_ORFs_that_Detected_a_Gene,Percentage_of_ORFs_that_Detected_a_Gene,Number_of_Genes_Detected,Percentage_of_Genes_Detected,Median_Length_of_All_ORFs,Median_Length_Difference,Minimum_Length_of_All_ORFs,Minimum_Length_Difference,Maximum_Length_of_All_ORFs,Maximum_Length_Difference,Median_GC_content_of_All_ORFs,Percent_Difference_of_All_ORFs_Median_GC,Median_GC_content_of_Matched_ORFs,Percent_Difference_of_Matched_ORF_GC,Number_of_ORFs_which_Overlap_Another_ORF,Percent_Difference_of_Overlapping_ORFs,Maximum_ORF_Overlap,Median_ORF_Overlap,Number_of_Matched_ORFs_Overlapping_Another_ORF,Percentage_Difference_of_Matched_Overlapping_CDSs,Maximum_Matched_ORF_Overlap,Median_Matched_ORF_Overlap,Number_of_Short-ORFs,Percent_Difference_of_Short-ORFs,Number_of_Short-Matched-ORFs,Percent_Difference_of_Short-Matched-ORFs,Number_of_Perfect_Matches,Percentage_of_Perfect_Matches,Number_of_Perfect_Starts,Percentage_of_Perfect_Starts,Number_of_Perfect_Stops,Percentage_of_Perfect_Stops,Number_of_Out_of_Frame_ORFs,Number_of_Matched_ORFs_Extending_a_Coding_Region,Percentage_of_Matched_ORFs_Extending_a_Coding_Region,Number_of_Matched_ORFs_Extending_Start_Region,Percentage_of_Matched_ORFs_Extending_Start_Region,Number_of_Matched_ORFs_Extending_Stop_Region,Percentage_of_Matched_ORFs_Extending_Stop_Region,Number_of_All_ORFs_on_Positive_Strand,Percentage_of_All_ORFs_on_Positive_Strand,Number_of_All_ORFs_on_Negative_Strand,Percentage_of_All_ORFs_on_Negative_Strand,Median_Start_Difference_of_Matched_ORFs,Median_Stop_Difference_of_Matched_ORFs,ATG_Start_Percentage,GTG_Start_Percentage,TTG_Start_Percentage,ATT_Start_Percentage,CTG_Start_Percentage,Other_Start_Codon_Percentage,TAG_Stop_Percentage,TAA_Stop_Percentage,TGA_Stop_Percentage,Other_Stop_Codon_Percentage,True_Positive,False_Positive,False_Negative,Precision,Recall,False_Discovery_Rate,Nucleotide_True_Positive,Nucleotide_False_Positive,Nucleotide_True_Negative,Nucleotide_False_Negative,Nucleotide_Precision,Nucleotide_Recall,Nucleotide_False_Discovery_Rate,ORF_Nucleotide_Coverage_of_Genome,Matched_ORF_Nucleotide_Coverage_of_Genome
995,109.03,190,19.10,190,39.92,335.0,-62.17,89,-21.24,3152,-41.81,31.50,0.20,32.83,4.42,279,26.24,135,0.00,36,-83.71,31,4.50,443,1826.09,19,-17.39,128,67.37,162,85.26,154,81.05,0,0,0.00,4,2.11,0,0.00,570,0.57,425,0.43,67.5,-85.5,63.12,15.28,21.61,0.00,0.00,0.00,11.06,27.44,61.51,0.00,0.40,1.69,0.60,0.19,0.40,0.81,0.82,0.31,0.69,0.18,0.96,0.82,0.04,77.15,24.47
CDS_Gene_Coverage_of_Genome:
90.62
Start_Position_Difference:
-78,33,93,294,144,408,3,18,156,-42,45,90,333,333,-39,111,201,93,120,-354,-150,-366,117,-138,-240,123,-153,-51
Stop_Position_Difference:
-192,-147,108,-216,87,-678,-96,-156,-321,-240,-168,-162,-51,-126,-33,-3,-93,-12,-204,-189,-156,237,-45,-219,-201,-537,-30,-78,159,243,60,21,15,183,288,6
Alternative_Starts_Predicted:

Alternative_Stops_Predicted:

Undetected_Gene_Metrics:
ATG_Start ,GTG_Start ,TTG_Start ,ATT_Start ,CTG_Start ,Alternative_Start_Codon ,TGA_Stop ,TAA_Stop ,TAG_Stop ,Alternative_Stop_Codon ,Median_Length ,ORFs_on_Positive_Strand ,ORFs_on_Negative_Strand
88.46,7.69,3.85,0.00,0.00,0.00,0.00,74.13,25.87,0.00,1047.50,156,130
Perfect_Match_Genes:
>Myco_686_1828_+ 
ATGAAAATATTAATTAATAAAAGTGAATTGAATAAAATTTTGAAAAAAATGAATAACGTTATTATTTCCAATAACAAAATAAAACCACATCATTCATATTTTTTAATAGAGGCAAAAGAAAAAGAAATAAACTTTTATGCTAACAATGAATACTTTTCTGTCAAATGTAATTTAAATAAAAATATTGATATTCTTGAACAAGGCTCCTTAATTGTTAAAGGAAAAATTTTTAACGATCTTATTAATGGCATAAAAGAAGAGATTATTACTATTCAAGAAAAAGATCAAACACTTTTGGTTAAAACAAAAAAAACAAGTATTAATTTAAACACAATTAATGTGAATGAATTTCCAAGAATAAGGTTTAATGAAAAAAACGATTTAAGTGAATTTAATCAATTCAAAATAAATTATTCACTTTTAGTAAAAGGCATTAAAAAAATTTTTCACTCAGTTTCAAATAATCGTGAAATATCTTCTAAATTTAATGGAGTAAATTTCAATGGATCCAATGGAAAAGAAATATTTTTAGAAGCTTCTGACACTTATAAACTATCTGTTTTTGAGATAAAGCAAGAAACAGAACCATTTGATTTCATTTTGGAGAGTAATTTACTTAGTTTCATTAATTCTTTTAATCCTGAAGAAGATAAATCTATTGTTTTTTATTACAGAAAAGATAATAAAGATAGCTTTAGTACAGAAATGTTGATTTCAATGGATAACTTTATGATTAGTTACACATCGGTTAATGAAAAATTTCCAGAGGTAAACTACTTTTTTGAATTTGAACCTGAAACTAAAATAGTTGTTCAAAAAAATGAATTAAAAGATGCACTTCAAAGAATTCAAACTTTGGCTCAAAATGAAAGAACTTTTTTATGCGATATGCAAATTAACAGTTCTGAATTAAAAATAAGAGCTATTGTTAATAATATCGGAAATTCTCTTGAGGAAATTTCTTGTCTTAAATTTGAAGGTTATAAACTTAATATTTCTTTTAACCCAAGTTCTCTATTAGATCACATAGAGTCTTTTGAATCAAATGAAATAAATTTTGATTTCCAAGGAAATAGTAAGTATTTTTTGATAACCTCTAAAAGTGAACCTGAACTTAAGCAAATATTGGTTCCTTCAAGATAA 

>Myco_4812_7322_+ 
ATGGCAAAGCAACAAGATCAAGTAGATAAGATTCGTGAAAACTTAGACAATTCAACTGTCAAAAGTATTTCATTAGCAAATGAACTTGAGCGTTCATTCATGGAATATGCTATGTCAGTTATTGTTGCTCGTGCTTTACCTGATGCTAGAGATGGACTTAAACCAGTTCATCGTCGTGTTCTTTATGGTGCTTATATTGGTGGCATGCACCATGATCGTCCTTTTAAAAAGTCTGCGAGGATTGTTGGTGATGTAATGAGTAAATTCCACCCTCATGGTGATATGGCAATATATGACACCATGTCAAGAATGGCTCAAGACTTTTCATTAAGATACCTTTTAATTGATGGTCATGGTAATTTTGGTTCTATAGATGGTGATAGACCTGCTGCACAACGTTATACAGAAGCAAGATTATCTAAACTTGCAGCAGAACTTTTAAAAGATATTGATAAAGATACAGTTGACTTTATTGCTAATTATGATGGTGAGGAAAAAGAACCAACTGTTCTACCAGCAGCTTTCCCTAACTTACTTGCAAATGGTTCTAGTGGGATTGCAGTTGGAATGTCAACATCTATTCCTTCCCATAATCTCTCTGAATTAATTGCGGGTTTAATCATGTTAATTGATAATCCTCAATGCACTTTTCAAGAATTATTAACTGTAATTAAAGGACCTGATTTTCCAACAGGAGCTAACATTATCTACACAAAAGGAATTGAAAGCTACTTTGAAACAGGTAAAGGCAATGTAGTAATTCGTTCTAAAGTTGAGATAGAACAATTGCAAACAAGAAGTGCATTAGTTGTAACTGAAATTCCTTACATGGTTAACAAAACTACCTTAATTGAAAAGATTGTAGAACTTGTTAAAGCTGAAGAGATTTCAGGAATTGCTGATATCCGTGATGAATCCTCTCGAGAAGGAATAAGGTTAGTGATTGAAGTAAAACGCGACACTGTACCTGAAGTTTTATTAAATCAACTTTTTAAATCAACAAGATTACAAGTACGCTTCCCTGTTAATATGCTTGCTTTAGTTAAAGGAGCTCCTGTACTTCTCAACATGAAACAAGCTTTGGAAGTATATCTTGATCATCAAATTGATGTTCTTGTTAGAAAAACAAAGTTTGTGCTTAATAAACAACAAGAACGTTATCACATTTTAAGCGGACTTTTAATTGCTGCTTTAAATATTGATGAGGTTGTTGCAATTATTAAAAAATCAGCAAATAACCAGGAAGCAATTAATACATTAAATACAAAGTTTAAGCTTGATGAAATTCAAGCTAAAGCAGTTCTTGACATGCGTTTAAGGAGCTTAAGCGTACTTGAAGTTAACAAACTTCAAACTGAACAAAAAGAGTTAAAAGATTCAATTGAATTTTGTAAGAAAGTGTTAGCTGATCAAAAATTACAGCTAAAAATAATCAAAGAGGAATTGCAAAAAATCAATGATCAGTTTGGTGATGAAAGAAGAAGTGAAATTCTCTATGATATCTCTGAGGAAATTGATGATGAATCATTGATAAAAGTTGAGAATGTAGTGATAACTATGTCTACAAATGGTTATCTAAAAAGGATTGGAGTTGATGCTTATAATCTTCAACATCGTGGTGGAGTTGGGGTTAAAGGGCTAACTACTTATGTTGATGATAGTATTAGTCAATTATTGGTCTGTTCAACTCACTCTGACTTATTATTTTTTACTGATAAGGGTAAGGTTTATAGAATTAGAGCTCATCAAATTCCCTATGGTTTTAGAACAAATAAAGGTATTCCCGCTGTTAACTTAATCAAAATTGAAAAGGATGAAAGAATTTGTTCATTGTTATCTGTTAATAACTATGATGATGGTTATTTCTTTTTCTGTACTAAAAATGGAATTGTTAAAAGAACGAGCTTGAATGAATTCATCAACATCTTAAGTAATGGTAAGCGGGCTATATCTTTTGATGATAATGACACTTTGTATTCAGTAATTAAAACCCACGGAAATGATGAGATTTTTATTGGTTCTACCAATGGATTTGTTGTTCGCTTCCATGAAAATCAACTCAGAGTTCTTTCAAGAACAGCAAGAGGTGTATTTGGTATCAGTTTAAATAAAGGAGAATTTGTTAATGGACTATCAACTTCAAGCAACGGTAGCTTACTTTTATCAGTCGGTCAAAATGGAATAGGTAAATTAACGAGCATAGATAAATATAGACTCACAAAACGTAATGCTAAGGGAGTTAAAACTCTAAGGGTTACTGATAGAACAGGCCCTGTTGTTACAACAACCACTGTTTTTGGTAATGAGGATCTTTTAATGATTTCCTCTGCTGGTAAAATTGTGCGTACCAGTTTACAAGAACTTTCAGAACAAGGTAAAAACACTTCTGGTGTTAAGTTAATTAGATTAAAAGATAATGAACGTTTAGAAAGAGTAACTATCTTTAAAGAAGAGTTAGAAGACAAAGAAATGCAACTAGAAGATGTTGGATCCAAACAAATTACGCAATAA 
.........
Partial_Match_Genes:
Gene:9923_11251_+_ATG_TAA 
ATGAAAAGCGAAATTAATATTTTTGCACTAGCAACTGCACCTTTTAATAGTGCATTACATATTATTAGGTTTTCTGGTCCTGATGTTTATGAGATTTTAAACAAGATAACTAATAAAAAAATAACAAGAAAAGGGATGCAAATTCAACGCACATGGATAGTTGATGAAAACAATAAGCGAATTGATGATGTGCTATTATTTAAATTTGTCTCTCCAAATTCTTATACAGGAGAAGATTTAATTGAAATTTCTTGTCATGGTAACATGTTGATCGTTAATGAAATTTGCGCACTTCTTTTAAAAAAAGGAGGTGTTTATGCCAAACCTGGTGAATTTACCCAAAGGAGTTTTTTAAATGGAAAAATGAGTTTACAACAAGCTAGTGCTGTAAATAAATTGATTTTATCTCCTAACTTATTAGTTAAAGATATAGTCTTAAATAATTTAGCGGGTGAAATGGATCAACAATTAGAACAAATAGCTCAACAAGTTAATCAATTAGTAATGCAAATGGAAGTAAACATTGATTATCCAGAATATCTTGATGAACAAGTAGAACTATCAACTTTAAATAATAAAGTTAAATTGATTATTGAAAAGCTTAAAAGAATTATTGAAAATAGTAAACAACTCAAAAAACTTCACGATCCTTTTAAAATTGCCATTATAGGCGAAACTAATGTAGGTAAATCTTCTTTACTCAACGCTTTATTAAATCAAGATAAAGCGATAGTTTCAAATATTAAAGGTAGTACACGCGATGTTGTTGAAGGGGATTTCAATTTAAATGGTTATTTAATCAAGATCTTAGATACTGCAGGTATCCGTAAACATAAAAGTGGGCTTGAAAAAGCAGGAATTAAAAAAAGCTTTGAATCTATAAAGCAAGCTAATTTGGTTATTTATCTTTTAGATGCAACACATCCAAAGAAAGATCTTGAATTAATTAGTTTTTTTAAGAAAAATAAAAAGGATTTTTTTGTTTTCTATAACAAAAAAGATTTAATTACAAATAAGTTTGAAAATAGTATTTCTGCAAAGCAAAAAGATATTAAAGAATTAGTTGATTTATTAACTAAATATATTAACGAGTTTTATAAAAAAATAGATCAAAAAATCTATCTGATTGAAAATTGACAGCAAATTTTAATTGAAAAAATTAAAGAACAATTAGAACAGTTTTTAAAGCAACAAAAAAAATATTTATTTTTCGATGTTTTAGTTACCCATCTAAGAGAAGCTCAACAAGATATTCTTAAACTACTAGGTAAGGATGTAGGTTTTGATTTAGTTAATGAAATTTTTAATAATTTTTGTTTAGGAAAATAA 
ORF:9923_11059_+_ATG_TGA 
ATGAAAAGCGAAATTAATATTTTTGCACTAGCAACTGCACCTTTTAATAGTGCATTACATATTATTAGGTTTTCTGGTCCTGATGTTTATGAGATTTTAAACAAGATAACTAATAAAAAAATAACAAGAAAAGGGATGCAAATTCAACGCACATGGATAGTTGATGAAAACAATAAGCGAATTGATGATGTGCTATTATTTAAATTTGTCTCTCCAAATTCTTATACAGGAGAAGATTTAATTGAAATTTCTTGTCATGGTAACATGTTGATCGTTAATGAAATTTGCGCACTTCTTTTAAAAAAAGGAGGTGTTTATGCCAAACCTGGTGAATTTACCCAAAGGAGTTTTTTAAATGGAAAAATGAGTTTACAACAAGCTAGTGCTGTAAATAAATTGATTTTATCTCCTAACTTATTAGTTAAAGATATAGTCTTAAATAATTTAGCGGGTGAAATGGATCAACAATTAGAACAAATAGCTCAACAAGTTAATCAATTAGTAATGCAAATGGAAGTAAACATTGATTATCCAGAATATCTTGATGAACAAGTAGAACTATCAACTTTAAATAATAAAGTTAAATTGATTATTGAAAAGCTTAAAAGAATTATTGAAAATAGTAAACAACTCAAAAAACTTCACGATCCTTTTAAAATTGCCATTATAGGCGAAACTAATGTAGGTAAATCTTCTTTACTCAACGCTTTATTAAATCAAGATAAAGCGATAGTTTCAAATATTAAAGGTAGTACACGCGATGTTGTTGAAGGGGATTTCAATTTAAATGGTTATTTAATCAAGATCTTAGATACTGCAGGTATCCGTAAACATAAAAGTGGGCTTGAAAAAGCAGGAATTAAAAAAAGCTTTGAATCTATAAAGCAAGCTAATTTGGTTATTTATCTTTTAGATGCAACACATCCAAAGAAAGATCTTGAATTAATTAGTTTTTTTAAGAAAAATAAAAAGGATTTTTTTGTTTTCTATAACAAAAAAGATTTAATTACAAATAAGTTTGAAAATAGTATTTCTGCAAAGCAAAAAGATATTAAAGAATTAGTTGATTTATTAACTAAATATATTAACGAGTTTTATAAAAAAATAGATCAAAAAATCTATCTGATTGAAAATTGA 

Gene:11251_12039_+_ATG_TAA 
ATGGAATACTTTGATGCACATTGTCATTTAAATTGTGAACCTTTACTGAGTGAAATTGAAAAAAGCATCGCTAATTTCAAATTAATTAATTTAAAAGCAAATGTTGTAGGTACAGATTTGGATAATTCTAAAATTGCTGTTGAATTAGCTAAAAAATATCCTGATCTTTTAAAAGCAACCATAGGTATCCATCCAAATGATGTTCATTTAGTTGATTTTAAAAAGACAAAAAAACAACTTAATGAACTATTAATAAATAACAGAAATTTCATAAGTTGTATTGGTGAATATGGTTTTGATTATCACTACACAACAGAATTTATTGAATTGCAAAACAAATTCTTTGAGATGCAATTTGAAATAGCTGAAACTAATAAATTGGTTCACATGCTTCATATTCGTGATGCTCATGAAAAAATTTATGAAATATTAACAAGATTAAAGCCAACTCAACCTGTGATTTTTCATTGTTTCAGTCAAGATATAAATATTGCTAAAAAGCTACTATCATTAAAAGATTTAAATATTGACATCTTCTTTTCTATCCCAGGGATAGTTACTTTTAAGAATGCTCAAGCATTACATGAAGCTTTAAAGATTATTCCTAGTGAATTACTTTTAAGTGAAACTGACTCACCGTGATTAACCCCTTCTCCTTTTCGAGGCAAAGTTAACTGACCTGAATATGTAGTTCATACTGTTAGCACTGTTGCTGAAATAAAAAAAATAGAAATTGCTGAAATGAAGCGAATTATTGTTAAAAATGCAAAAAAATTATTTTGACATTAA 
ORF:11251_11892_+_ATG_TGA 
ATGGAATACTTTGATGCACATTGTCATTTAAATTGTGAACCTTTACTGAGTGAAATTGAAAAAAGCATCGCTAATTTCAAATTAATTAATTTAAAAGCAAATGTTGTAGGTACAGATTTGGATAATTCTAAAATTGCTGTTGAATTAGCTAAAAAATATCCTGATCTTTTAAAAGCAACCATAGGTATCCATCCAAATGATGTTCATTTAGTTGATTTTAAAAAGACAAAAAAACAACTTAATGAACTATTAATAAATAACAGAAATTTCATAAGTTGTATTGGTGAATATGGTTTTGATTATCACTACACAACAGAATTTATTGAATTGCAAAACAAATTCTTTGAGATGCAATTTGAAATAGCTGAAACTAATAAATTGGTTCACATGCTTCATATTCGTGATGCTCATGAAAAAATTTATGAAATATTAACAAGATTAAAGCCAACTCAACCTGTGATTTTTCATTGTTTCAGTCAAGATATAAATATTGCTAAAAAGCTACTATCATTAAAAGATTTAAATATTGACATCTTCTTTTCTATCCCAGGGATAGTTACTTTTAAGAATGCTCAAGCATTACATGAAGCTTTAAAGATTATTCCTAGTGAATTACTTTTAAGTGAAACTGACTCACCGTGA 
.......
Missed_Genes:
>Myco_1828_2760_+ 
ATGAATCTTTACGATCTTTTAGAACTACCAACTACAGCATCAATAAAAGAAATAAAAATTGCTTATAAAAGATTAGCAAAGCGTTATCACCCTGATGTAAATAAATTAGGTTCGCAAACTTTTGTTGAAATTAATAATGCTTATTCAATATTAAGTGATCCTAACCAAAAGGAAAAATATGATTCAATGCTGAAAGTTAATGATTTTCAAAATCGCATCAAAAATTTAGATATTAGTGTTAGATGACATGAAAATTTCATGGAAGAACTCGAACTTCGTAAGAACTGAGAATTTGATTTTTTTTCATCTGATGAAGATTTCTTTTATTCTCCATTTACAAAAAACAAATATGCTTCCTTTTTAGATAAAGATGTTTCTTTAGCTTTTTTTCAGCTTTACAGCAAGGGCAAAATAGATCATCAATTGGAAAAATCTTTATTGAAAAGAAGAGATGTAAAAGAAGCTTGTCAACAGAATAAAAATTTTATTGAAGTTATAAAAGAGCAATATAACTATTTTGGTTGAATTGAAGCTAAGCGTTATTTCAATATTAATGTTGAACTTGAGCTCACACAGAGAGAGATAAGAGATAGAGATGTTGTTAACCTACCTTTAAAAATTAAAGTTATTAATAATGATTTTCCAAATCAACTCTGATATGAAATTTATAAAAACTATTCATTTCGCTTATCTTGAGATATAAAAAATGGTGAAATTGCTGAATTTTTCAATAAAGGTAATAGAGCTTTAGGATGAAAAGGTGACTTAATTGTCAGAATGAAAGTAGTTAATAAAGTAAACAAAAGACTGCGTATTTTTTCAAGCTTTTTTGAGAACGATAAATCTAAATTATGGTTCCTTGTTCCAAACGATAAACAAAGTAATCCTAATAAGGGCGTTTTTAACTATAAAACTCAGCACTTTATTGATTAA 

>Myco_2845_4797_+ 
ATGGAAGAAAATAACAAAGCAAATATCTATGACTCTAGTAGCATTAAGGTCCTTGAAGGACTTGAGGCTGTTAGAAAACGCCCTGGAATGTACATTGGTTCTACTGGCGAAGAAGGTTTGCATCACATGATCTGAGAGATAGTAGACAACTCAATTGATGAAGCAATGGGAGGTTTTGCCAGTTTTGTTAAGCTTACCCTTGAAGATAATTTTGTTACCCGTGTAGAGGATGATGGAAGAGGGATACCTGTTGATATCCATCCTAAGACTAATCGTTCTACAGTTGAAACAGTTTTTACAGTTCTACACGCTGGCGGTAAATTTGATAACGATAGCTATAAAGTGTCAGGTGGTTTACACGGTGTTGGTGCATCAGTTGTTAATGCGCTTAGTTCTTCTTTTAAAGTTTGAGTTTTTCGTCAAAATAAAAAGTATTTTCTCAGCTTTAGCGATGGAGGAAAGGTAATTGGAGATTTGGTCCAAGAAGGTAACTCTGAAAAAGAGCATGGAACAATTGTTGAGTTTGTTCCTGATTTCTCTGTAATGGAAAAGAGTGATTACAAACAAACTGTAATTGTAAGCAGACTCCAGCAATTAGCTTTTTTAAACAAGGGAATAAGAATTGACTTTGTTGATAATCGTAAACAAAACCCACAGTCTTTTTCTTGAAAATATGATGGGGGATTGGTTGAATATATCCACCACCTAAACAACGAAAAAGAACCACTTTTTAATGAAGTTATTGCTGATGAAAAAACTGAAACTGTAAAAGCTGTTAATCGTGATGAAAACTACACAGTAAAGGTTGAAGTTGCTTTTCAATATAACAAAACATACAACCAATCAATTTTCAGTTTTTGTAACAACATTAATACTACAGAAGGTGGAACCCATGTGGAAGGTTTTCGTAATGCACTTGTTAAGATCATTAATCGCTTTGCTGTTGAAAATAAATTCCTAAAAGATAGTGATGAAAAGATTAACCGTGATGATGTTTGTGAAGGATTAACTGCTATTATTTCCATTAAACACCCAAACCCACAATATGAAGGACAAACTAAAAAGAAGTTAGGTAATACTGAGGTAAGACCTTTAGTTAATAGTGTTGTTAGTGAAATCTTTGAACGCTTCATGTTAGAAAACCCACAAGAAGCAAACGCTATCATCAGAAAAACACTTTTAGCTCAAGAAGCGAGAAGAAGAAGTCAAGAGGCTAGGGAGTTAACTCGTCGTAAATCACCTTTTGATAGTGGTTCATTACCAGGTAAATTAGCTGATTGTACAACCAGAGATCCTTCGATTAGTGAACTTTACATTGTTGAGGGTGATAGTGCTGGTGGCACTGCTAAAACAGGAAGAGATCGTTATTTTCAAGCTATCTTACCCTTAAGAGGAAAGATTTTAAACGTTGAAAAATCTAACTTTGAACAAATCTTTAATAATGCAGAAATTTCTGCATTAGTGATGGCAATAGGCTGTGGGATTAAACCTGATTTTGAACTTGAAAAACTTAGATATAGCAAGATTGTGATCATGACAGATGCTGATGTTGATGGTGCACACATAAGAACACTTCTCTTAACTTTCTTTTTTCGCTTTATGTATCCTTTGGTTGAACAAGGCAATATTTTTATTGCTCAACCCCCACTTTATAAAGTGTCATATTCCCATAAGGATTTATACATGCACACTGATGTTCAACTTGAACAGTGAAAAAGTCAAAACCCTAACGTAAAGTTTGGGTTACAAAGATATAAAGGACTTGGAGAAATGGATGCATTGCAGCTGTGAGAAACAACAATGGATCCTAAGGTTAGAACATTGTTAAAAGTTACTGTTGAAGATGCTTCTATTGCTGATAAAGCTTTTTCACTGTTGATGGGTGATGAAGTTCCCCCAAGAAGAGAATTTATTGAAAAAAATGCTCGTAGTGTTAAAAACATTGATATTTAA 

>Myco_7294_8547_+ 
ATGTTGGATCCAAACAAATTACGCAATAACTATGATTTCTTTAAAAAGAAACTGTTAGAAAGAAATGTAAATGAGCAATTATTAAATCAGTTTATTCAAACTGATAAACTAATGCGCAAAAACTTGCAACAACTTGAACTTGCTAACCAAAAACAAAGCTTGTTGGCAAAACAAGTTGCTAAGCAAAAAGATAATAAAAAGCTATTAGCTGAATCAAAAGAACTTAAGCAGAAGATTGAAAACTTAAATAATGCTTATAAAGATTCACAAAACATTAGTCAAGATTTACTTCTAAATTTTCCTAATATTGCTCATGAATCAGTTCCTGTTGGTAAAAATGAATCAGCAAACTTAGAACTTCTTAAAGAAGGGAGAAAACCAGTTTTTGATTTCAAACCTTTACCACATCGAGAGTTATGTGAAAAGTTAAATTTAGTTGCTTTTGATAAAGCTACTAAGATTAGTGGAACTAGGTTTGTTGCATATACAGATAAAGCAGCTAAACTACTTAGAGCGATAACTAATCTAATGATTGACCTTAATAAAAGCAAGTATCAAGAATGAAACCTGCCAGTTGTTATTAATGAATTAAGTTTAAGATCAACCGGACAACTACCTAAGTTTAAAGATGATGTTTTTAAACTAGAAAACACCCGTTATTATCTTTCTCCAACTTTAGAGGTACAACTTATCAATTTACATGCTAATGAAATTTTTAATGAAGAAGATTTACCTAAATACTACACTGCAACAGGTATTAACTTTCGTCAAGAAGCGGGTAGTGCTGGTAAACAAACCAAAGGAACTATTAGATTGCATCAGTTTCAAAAAACTGAGTTAGTTAAGTTTTGTAAACCTGAAAATGCTATCAATGAATTGGAAGCAATGGTTAGAGATGCTGAACAAATCTTAAAGGCACTTAAGTTACCTTTTAGAAGGTTATTGTTATGTACTGGTGATATGGGCTTTAGTGCTGAAAAAACATATGATCTTGAAGTTTGAATGGCAGCTAGCAATGAATATCGTGAAGTTTCTTCTTGTTCATCTTGTGGTGATTTTCAAGCAAGAAGAGCTATGATTCGTTACAAAGATATTAACAACGGTAAAAACAGTTATGTTGCTACTTTAAATGGAACAGCATTATCTATTGATAGAATTTTTGCTGCAATTCTAGAAAATTTTCAAACAAAAGATGGCAAAATTCTTATCCCACAAGCATTAAAAAAATACCTTGATTTTGACACAATCAAGTAA 
......
 
ORFs_Without_Corresponding_Gene_In_Reference_Metrics:
ATG_Start ,GTG_Start ,TTG_Start ,ATT_Start ,CTG_Start ,Alternative_Start_Codon ,TGA_Stop ,TAA_Stop ,TAG_Stop ,Alternative_Stop_Codon ,Median_Length ,ORFs_on_Positive_Strand ,ORFs_on_Negative_Strand
58.39,17.14,24.47,0.00,0.00,0.00,71.55,20.62,7.83,0.00,287.00,449,356
ORF_Without_Corresponding_Gene_in_Reference:
>Prodigal_1828_2073_+ 
ATGAATCTTTACGATCTTTTAGAACTACCAACTACAGCATCAATAAAAGAAATAAAAATTGCTTATAAAAGATTAGCAAAGCGTTATCACCCTGATGTAAATAAATTAGGTTCGCAAACTTTTGTTGAAATTAATAATGCTTATTCAATATTAAGTGATCCTAACCAAAAGGAAAAATATGATTCAATGCTGAAAGTTAATGATTTTCAAAATCGCATCAAAAATTTAGATATTAGTGTTAGATGA
>Prodigal_2605_2760_+ 
ATGAAAGTAGTTAATAAAGTAAACAAAAGACTGCGTATTTTTTCAAGCTTTTTTGAGAACGATAAATCTAAATTATGGTTCCTTGTTCCAAACGATAAACAAAGTAATCCTAATAAGGGCGTTTTTAACTATAAAACTCAGCACTTTATTGATTAA
>Prodigal_2845_2979_+ 
ATGGAAGAAAATAACAAAGCAAATATCTATGACTCTAGTAGCATTAAGGTCCTTGAAGGACTTGAGGCTGTTAGAAAACGCCCTGGAATGTACATTGGTTCTACTGGCGAAGAAGGTTTGCATCACATGATCTGA
>Prodigal_3010_3255_+ 
ATGGGAGGTTTTGCCAGTTTTGTTAAGCTTACCCTTGAAGATAATTTTGTTACCCGTGTAGAGGATGATGGAAGAGGGATACCTGTTGATATCCATCCTAAGACTAATCGTTCTACAGTTGAAACAGTTTTTACAGTTCTACACGCTGGCGGTAAATTTGATAACGATAGCTATAAAGTGTCAGGTGGTTTACACGGTGTTGGTGCATCAGTTGTTAATGCGCTTAGTTCTTCTTTTAAAGTTTGA
>Prodigal_3319_3513_+ 
TTGGTCCAAGAAGGTAACTCTGAAAAAGAGCATGGAACAATTGTTGAGTTTGTTCCTGATTTCTCTGTAATGGAAAAGAGTGATTACAAACAAACTGTAATTGTAAGCAGACTCCAGCAATTAGCTTTTTTAAACAAGGGAATAAGAATTGACTTTGTTGATAATCGTAAACAAAACCCACAGTCTTTTTCTTGA
>Prodigal_3529_4557_+ 
TTGGTTGAATATATCCACCACCTAAACAACGAAAAAGAACCACTTTTTAATGAAGTTATTGCTGATGAAAAAACTGAAACTGTAAAAGCTGTTAATCGTGATGAAAACTACACAGTAAAGGTTGAAGTTGCTTTTCAATATAACAAAACATACAACCAATCAATTTTCAGTTTTTGTAACAACATTAATACTACAGAAGGTGGAACCCATGTGGAAGGTTTTCGTAATGCACTTGTTAAGATCATTAATCGCTTTGCTGTTGAAAATAAATTCCTAAAAGATAGTGATGAAAAGATTAACCGTGATGATGTTTGTGAAGGATTAACTGCTATTATTTCCATTAAACACCCAAACCCACAATATGAAGGACAAACTAAAAAGAAGTTAGGTAATACTGAGGTAAGACCTTTAGTTAATAGTGTTGTTAGTGAAATCTTTGAACGCTTCATGTTAGAAAACCCACAAGAAGCAAACGCTATCATCAGAAAAACACTTTTAGCTCAAGAAGCGAGAAGAAGAAGTCAAGAGGCTAGGGAGTTAACTCGTCGTAAATCACCTTTTGATAGTGGTTCATTACCAGGTAAATTAGCTGATTGTACAACCAGAGATCCTTCGATTAGTGAACTTTACATTGTTGAGGGTGATAGTGCTGGTGGCACTGCTAAAACAGGAAGAGATCGTTATTTTCAAGCTATCTTACCCTTAAGAGGAAAGATTTTAAACGTTGAAAAATCTAACTTTGAACAAATCTTTAATAATGCAGAAATTTCTGCATTAGTGATGGCAATAGGCTGTGGGATTAAACCTGATTTTGAACTTGAAAAACTTAGATATAGCAAGATTGTGATCATGACAGATGCTGATGTTGATGGTGCACACATAAGAACACTTCTCTTAACTTTCTTTTTTCGCTTTATGTATCCTTTGGTTGAACAAGGCAATATTTTTATTGCTCAACCCCCACTTTATAAAGTGTCATATTCCCATAAGGATTTATACATGCACACTGATGTTCAACTTGAACAGTGA
....
ORFs_Which_Detected_more_than_one_Gene:

```


## GFF Tools:

### GFF-Adder:

GFF-Adder allows for the addition of predicted CDSs to an existing reference annotation (GFF or another tool) which produces a new GFF containing the original
genes plus the new CDS from another prediction. Default filtering will remove additional CDSs that overlap existing genes by more than 50 nt.
The ```-gi``` option can be used to allow for different genomic elements to be accounted for, other than only CDSs in the reference annotation.

For Help: ```GFF-Adder -h ```

```python
Thank you for using ORForise
Please report any issues to: https://github.com/NickJD/ORForise/issues
#####
usage: GFF_Adder.py [-h] -dna GENOME_DNA -ref REFERENCE_ANNOTATION -at ADDITIONAL_TOOL -add ADDITIONAL_ANNOTATION -o
                    OUTPUT_FILE [-rt REFERENCE_TOOL] [-gi GENE_IDENT] [-gene_ident GENE_IDENT] [-olap OVERLAP]

ORForise v1.4.1: GFF-Adder Run Parameters.

Required Arguments:
  -dna GENOME_DNA       Genome DNA file (.fa) which both annotations are based on
  -ref REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -at ADDITIONAL_TOOL   Which format to use for additional annotation?
  -add ADDITIONAL_ANNOTATION
                        Which annotation file to add to reference annotation?
  -o OUTPUT_FILE        Output filename

Optional Arguments:
  -rt REFERENCE_TOOL    Which tool format to use as reference? - If not provided, will default to standard Ensembl
                        GFF format, can be Prodigal or any of the other tools available
  -gi GENE_IDENT        Identifier used for extraction of "genic" regions from reference annotation "CDS,rRNA,tRNA":
                        Default for is "CDS"
  -gene_ident GENE_IDENT
                        Identifier used for identifying genomic features in reference annotation "CDS,rRNA,tRNA"
  -olap OVERLAP         Maximum overlap between reference and additional genic regions (CDS,rRNA etc) - Default: 50
                        nt
```

#### Example: Running GFF-Adder to combine the additional CDS predictions made by Prodial to the canonical annotations from Ensembl.
``` GFF-Adder -dna ~/Testing/Myco.fa -ref ~/Testing/Myco.gff  -at Prodigal -add ~/Testing/Prodigal_Myco.gff -o ~/Testing/Myco_Ensembl_GFF_Adder_Prodigal.gff ```
#### Example Output: ~/ORForise/Testing/Myco_Ensembl_GFF_Adder_Prodigal.gff
```
##gff-version	3
#	GFF-Adder
#	Run Date:2021-11-10
##Genome DNA File:./Testing/Myco.fa
##Original File: ./Testing/Myco.gff
##Additional File: ./Testing/Prodigal_Myco.gff
.......
Chromosome	Reference_Annotation	CDS	68522	70225	.	-	.	ID=Original_Annotation
Chromosome	Reference_Annotation	CDS	70530	72572	.	+	.	ID=Original_Annotation
Chromosome	Reference_Annotation	CDS	72523	73434	.	+	.	ID=Original_Annotation
Chromosome	Prodigal	CDS	73445	73648	.	+	.	ID=Additional_Annotation
Chromosome	Reference_Annotation	CDS	73690	77685	.	+	.	ID=Original_Annotation
Chromosome	Reference_Annotation	CDS	77685	79085	.	+	.	ID=Original_Annotation
Chromosome	Reference_Annotation	CDS	79089	81035	.	+	.	ID=Original_Annotation
Chromosome	Reference_Annotation	CDS	81046	82596	.	+	.	ID=Original_Annotation
Chromosome	Reference_Annotation	CDS	82620	84044	.	+	.	ID=Original_Annotation
Chromosome	Prodigal	CDS	84082	84312	.	+	.	ID=Additional_Annotation
Chromosome	Prodigal	CDS	84532	84744	.	-	.	ID=Additional_Annotation
Chromosome	Prodigal	CDS	84776	85051	.	+	.	ID=Additional_Annotation
```

### GFF-Intersector:

GFF-Intersector enables the aggregation of different genome annotations and CDS predictions and creates a single GFF
representing the intersection of the two existing annotations.
GFF-Intersector also provides an option to allow the retention of genes that have a user defined difference (minimum % coverage and in-frame).
The ```-gi``` option can be used to allow for different genomic elements to be accounted for, other than only CDSs in the reference annotation.

For Help: ```GFF-Intersector -h ``` 
```python
Thank you for using ORForise
Please report any issues to: https://github.com/NickJD/ORForise/issues
#####
usage: GFF_Intersector.py [-h] -dna GENOME_DNA -ref REFERENCE_ANNOTATION -at ADDITIONAL_TOOL -add
                          ADDITIONAL_ANNOTATION -o OUTPUT_FILE [-rt REFERENCE_TOOL] [-gi GENE_IDENT] [-cov COVERAGE]

ORForise v1.4.1: GFF-Intersector Run Parameters.

Required Arguments:
  -dna GENOME_DNA       Genome DNA file (.fa) which both annotations are based on
  -ref REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -at ADDITIONAL_TOOL   Which format to use for additional annotation?
  -add ADDITIONAL_ANNOTATION
                        Which annotation file to add to reference annotation?
  -o OUTPUT_FILE        Output filename

Optional Arguments:
  -rt REFERENCE_TOOL    Which tool format to use as reference? - If not provided, will default to standard Ensembl
                        GFF format, can be Prodigal or any of the other tools available
  -gi GENE_IDENT        Identifier used for extraction of "genic" regions from reference annotation "CDS,rRNA,tRNA":
                        Default for is "CDS"
  -cov COVERAGE         Percentage coverage of reference annotation needed to confirm intersection - Default: 100 ==
                        exact match
```

#### Example: Running GFF-Intersector to combine the additional CDS predictions made by Prodial to the canonical annotations from Ensembl.
``` GFF-Intersector -dna ~/Testing/Myco.fa -ref ~/Testing/Myco.gff -at Prodigal -add ~/Testing/Prodigal_Myco.gff -o ~/Testing/Myco_Ensembl_GFF_Intersector_Prodigal.gff```

#### Example Output: ~/Testing/Myco_Ensembl_GFF_Intersector_Prodigal.gff
```
##gff-version	3
#	GFF-Intersector
#	Run Date:2021-11-10
##Genome DNA File:./Testing/Myco.fa
##Original File: ./Testing/Myco.gff
##Intersecting File: ./Testing/Prodigal_Myco.gff
Chromosome	original	CDS	686	1828	.	+	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	4812	7322	.	+	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	8551	9183	.	+	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	22389	23558	.	+	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	29552	30124	.	+	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	31705	32325	.	-	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	49376	49642	.	+	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	59082	59753	.	+	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	61014	61406	.	+	.	ID=Original_Annotation;Coverage=100
Chromosome	original	CDS	82620	84044	.	+	.	ID=Original_Annotation;Coverage=100
```


# Genomes Available:

The .fa and .gff files (from Ensembl Bacteria Release 46) below are available in the Genomes directory.

* *Bacillus subtilis* - Strain BEST7003 - Assembly ASM52304v1
* *Caulobacter crescentus* - Strain CB15 - Assembly ASM690v1
* *Escherichia coli K-12* - Strain ER3413 - Assembly ASM80076v1
* *Mycoplasma genitalium* - Strain G37 - Assembly ASM2732v1
* *Pseudomonas fluorescens* - Strain UK4 - Assembly ASM73042v1
* *Staphylococcus aureus* - Strain 502A - Assembly ASM59796v1

# Prediction Tools Available:

There are two Groups of tools - Those which do require a pre-built model and those which do not. \
For the example runs provided, each tool is listed with the non-default options used and their predictions for each of the 6 model organisms are available in their respective
directories.
ORForise only needs the tool name and the annotation file produced from any available model to undertake the analysis.

## GFF Standard Format:

The GFF Tool directory allows for the analysis of user-provided annotations in the standard GFF3 format. \
This can be used to compare different cannonical annotations with eachother or additional tools which use the GFF3 format.

## Model Based Tools:

**Augustus - Version 3.3.3** - http://bioinf.uni-greifswald.de/augustus/  
This tool has three comparisons with the organism models *E. coli* and *S. aureus* and *H. sapiens*.

**EasyGene - Version 1.2** - http://www.cbs.dtu.dk/services/EasyGene/  
This tool has two comparisons with the organism models *E. coli - K12* and *S. aureus Mu50*.

**FGENESB - Version '2020'** - http://www.softberry.com/berry.phtml?topic=fgenesb&group=programs&subgroup=gfindb  
This tool has two comparisons with the organism models *E. coli - K12* and *S. aureus MU50*.

**GeneMark - Version 2.5** - http://exon.gatech.edu/GeneMark/gm.cgi  
This tool has two comparisons with the organism models *E. coli - K12 - MG165* and *S. aureus Mu50*.

**GeneMark.hmm - Version 3.2.5** -  http://exon.gatech.edu/GeneMark/gmhmmp.cgi  
This tool has two comparisons with the organism models *E. coli - K12 - MG165* and *S. aureus Mu50*.

## Self-Training/Non-Model Based Tools

**FragGeneScan - Version 1.3.0** - https://omics.informatics.indiana.edu/FragGeneScan/    
The 'complete' genome option was selected and GFF was chosen as output type.

**GeneMark HA - Version 3.25** - http://exon.gatech.edu/GeneMark/heuristic_gmhmmp.cgi  
GFF was chosen as output type.

**GeneMarkS - Version 4.25** - http://exon.gatech.edu/GeneMark/genemarks.cgi  
GFF was chosen as output type.

**GeneMarkS-2 - Version '2020'** - http://exon.gatech.edu/GeneMark/genemarks2.cgi   
GFF3 was chosen as output type.

**GLIMMER-3 - Version 3.02** - http://ccb.jhu.edu/software/glimmer/index.shtml  
Default parameters from manual were used.

**MetaGene - Version 2.24.0** - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1636498/  
Default options were used.

**MetaGeneAnnotator - Version '2008/8/19'** - http://metagene.nig.ac.jp/  
Defaults options were used.

**MetaGeneMark - Version '2020'** - http://exon.gatech.edu/meta_gmhmmp.cgi  
GFF was chosen as output type.

**Prodigal - Version 2.6.3** - https://github.com/hyattpd/Prodigal  
GFF was chosen as output type.

**TransDecoder - Version 5.5.0** - https://github.com/TransDecoder/TransDecoder/wiki  
Defaults options were used.

**Balrog - Version 2021`** - https://github.com/salzberg-lab/Balrog  
Defaults options were used.


