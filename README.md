# ORForise - Prokaryote Genome Annotation Comparison and Analysis Platform

Platform for analysing Prokaryote CoDing Sequence (CDS) Gene Predictors. \
Novel genome annotations can be compared to a provided reference annotation from Ensembl (or any given GFF annotation) 
and predictions from other tools.

## Requirements and Installation:

The ORForise platform is written in Python3.8 and only requires the NumPy library which is standard in most base
installations of Python3. \
Usually, ```pip3 install numpy``` is adequate to install NumPy.

### Intallation:

The ORForise platform is available via github ```git clone https://github.com/NickJD/ORForise``` and the pip Python package manager ```pip3 install ORForise```. \
For both methods of 'installation', cloning the repository or installing ORForise via pip, the same input files are required when running the code as shown in the examples below. 

To run, you need:
* Input Genome FASTA and corresponding GFF file or CDS predictions with the annotated genes for the genome you want to use as reference.
* A prediction output from one of the compatible tools for the same genome.

### How to add your own Genome:

Corresponding FASTA and GFF files must be provided for the genome the analysis is to be performed on, including the corresponding output of any tools to compare.

### How to add your own tool:

If the new tool reports its predictions in GFF you can present ORForise with "GFF" for either the reference ```-rt``` or prediction ```-t``` option.
If the tool uses another non-standard format, a request can be made to add it as an option via GitHub.

## CDS Prediction Analysis:

### Use-cases: (Running if via pip)

For Help: ```python python3 -m ORForise.Annotation_Compare -h ```

```python
usage: Annotation_Compare.py [-h] -dna GENOME_DNA [-rt REFERENCE_TOOL] -ref REFERENCE_ANNOTATION -t TOOL -tp TOOL_PREDICTION [-o OUTNAME] [-v {True,False}]

optional arguments:
  -h, --help            show this help message and exit
  -dna GENOME_DNA, --genome_DNA GENOME_DNA
                        Genome DNA file (.fa) which both annotations are based on
  -rt REFERENCE_TOOL, --reference_tool REFERENCE_TOOL
                        What type of Annotation to compare to? -- Leave blank for Ensembl reference- Provide tool name to compare output from two tools (GeneMarkS)
  -ref REFERENCE_ANNOTATION, --reference_annotation REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -t TOOL, --tool TOOL  Which tool to analyse? (Prodigal)
  -tp TOOL_PREDICTION, --tool_prediction TOOL_PREDICTION
                        Tool genome prediction file (.gff) - Different Tool Parameters are compared individually via separate files
  -o OUTNAME, --outname OUTNAME
                        Define full output filename (format is CSV) - If not provided, summary will be printed to std-out
  -v {True,False}, --verbose {True,False}
                        Default - False: Print out runtime status
```

### Compare a novel genome annotation to an Ensembl annotation:

Genome annotation is a difficult process, even for Prokaryotes. ORForise allows the direct and systematic analysis of
a novel CDS prediction from a wide selection of tools to a reference Genome Annotation, such as those provided by
Ensembl Bacteria.

#### Example: Installation through pip will allow user to call the programs directly from the ORForise package.
```python
 python3 -m ORForise.Annotation_Compare -dna Genomes/Myco.fa -ref Genomes/Myco.gff -t Prodigal -tp Tools/Prodigal/Prodigal_Myco.gff
```
### Compare different novel annotations with each other on a single Genome:

If a reference Genome Annotation is not available or a direct comparison between two or more tools is wanted,
ORForise can be used as the example below.

## Aggregate CDS Prediction Analysis:

### Use-cases: (Running if via pip)

For Help: ```python python3 -m ORForise.Aggregate_Compare -h ```

```python

usage: Aggregate_Compare.py [-h] -dna GENOME_DNA -t TOOLS -tp TOOL_PREDICTIONS [-rt REFERENCE_TOOL] -ref REFERENCE_ANNOTATION [-o OUTNAME] [-v {True,False}]

optional arguments:
  -h, --help            show this help message and exit
  -dna GENOME_DNA, --genome_DNA GENOME_DNA
                        Genome DNA file (.fa) which both annotations are based on
  -t TOOLS, --tools TOOLS
                        Which tools to analyse? (Prodigal,GeneMarkS)
  -tp TOOL_PREDICTIONS, --tool_predictions TOOL_PREDICTIONS
                        Tool genome prediction file (.gff) - Providefile locations for each tool comma separated
  -rt REFERENCE_TOOL, --reference_tool REFERENCE_TOOL
                        What type of Annotation to compare to? -- Leave blank for Ensembl reference- Provide tool name to compare output from two tools
                        (GeneMarkS)
  -ref REFERENCE_ANNOTATION, --reference_annotation REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -o OUTNAME, --outname OUTNAME
                        Define full output filename (format is CSV) - If not provided, summary will be printed to std-out
  -v {True,False}, --verbose {True,False}
                        Default - False: Print out runtime status
```

#### Example: 
```python
python3 -m ORForise.Aggregate_Compare -ref /home/nick/Git/ORForise/src/Genomes/Myco.gff -dna /home/nick/Git/ORForise/src/Genomes/Myco.fa -t Prodigal,TransDecoder,GLIMMER_3 -tp /home/nick/Git/ORForise/src/ORForise/Tools/Prodigal/Prodigal_Myco.gff,/home/nick/Git/ORForise/src/ORForise/Tools/TransDecoder/TransDecoder_Myco.gff,/home/nick/Git/ORForise/src/ORForise/Tools/TransDecoder/TransDecoder_Myco.gff```
```
This will compare the Aggregate the predictions of Prodigal, TransDecoder and GLIMMER 3 against the Mycoplasma reference annotation provided by
Ensembl Bacteria.

## Annotation Comparison Output - The output format is the same for Annotation_Compare and Aggregate_Compare:
### Print to screen example - Prodigal prediction compared to Ensembl Bacteria reference annotation of *Escherichia coli*:
```bash
/usr/bin/python3.8 /home/nick/Git/ORForise/src/ORForise/Annotation_Compare.py -ref /home/nick/Git/ORForise/src/Genomes/E-coli.gff -dna /home/nick/Git/ORForise/src/Genomes/E-coli.fa -t Prodigal -tp /home/nick/Git/ORForise/src/ORForise/Tools/Prodigal/Prodigal_E-coli.gff -o /home/nick/Git/ORForise/src/ORForise/Tools/Prodigal/Prodigal_E-coli.csv
Genome Used: E-coli
Reference Used: /home/nick/Git/ORForise/src/Genomes/E-coli.gff
Tool Compared: Prodigal
Perfect Matches:3737[4052]
Partial Matches:236[4052]
Missed Genes:79[4052]
Complete
```
This is the default output of the comparison tools. 

### '-o' Example output to CSV file - Prodigal prediction compared to Ensembl Bacteria reference annotation of *Escherichia coli*:
The output is designed to be human-readable and interpretable by the included 'ORForise_Analysis' scripts. 
The example below presents the 12 'Representative' and 72 'All' Metrics but only shows one entry for each of the induvidual prediction reports (Perfect_Match_Genes,Partial_Match_Genes,Missed_Genes,Predicted_CDS_Without_Corresponding_Gene_in_Reference,Predicted_CDSs_Which_Detected_more_than_one_Gene).

```csv
Representative_Metrics:
Percentage_of_Genes_Detected,Percentage_of_ORFs_that_Detected_a_Gene,Percent_Difference_of_All_ORFs,Median_Length_Difference,Percentage_of_Perfect_Matches,Median_Start_Difference_of_Matched_ORFs,Median_Stop_Difference_of_Matched_ORFs,Percentage_Difference_of_Matched_Overlapping_CDSs,Percent_Difference_of_Short-Matched-ORFs,Precision,Recall,False_Discovery_Rate
98.05,93.15,5.21,-2.83,94.11,-3.0,N/A,-3.83,-16.63,0.93,0.98,0.07
All_Metrics:
Number_of_ORFs,Percent_Difference_of_All_ORFs,Number_of_ORFs_that_Detected_a_Gene,Percentage_of_ORFs_that_Detected_a_Gene,Number_of_Genes_Detected,Percentage_of_Genes_Detected,Median_Length_of_All_ORFs,Median_Length_Difference,Minimum_Length_of_All_ORFs,Minimum_Length_Difference,Maximum_Length_of_All_ORFs,Maximum_Length_Difference,Median_GC_content_of_All_ORFs,Percent_Difference_of_All_ORFs_Median_GC,Median_GC_content_of_Matched_ORFs,Percent_Difference_of_Matched_ORF_GC,Number_of_ORFs_which_Overlap_Another_ORF,Percent_Difference_of_Overlapping_ORFs,Maximum_ORF_Overlap,Median_ORF_Overlap,Number_of_Matched_ORFs_Overlapping_Another_ORF,Percentage_Difference_of_Matched_Overlapping_CDSs,Maximum_Matched_ORF_Overlap,Median_Matched_ORF_Overlap,Number_of_Short-ORFs,Percent_Difference_of_Short-ORFs,Number_of_Short-Matched-ORFs,Percent_Difference_of_Short-Matched-ORFs,Number_of_Perfect_Matches,Percentage_of_Perfect_Matches,Number_of_Perfect_Starts,Percentage_of_Perfect_Starts,Number_of_Perfect_Stops,Percentage_of_Perfect_Stops,Number_of_Out_of_Frame_ORFs,Number_of_Matched_ORFs_Extending_a_Coding_Region,Percentage_of_Matched_ORFs_Extending_a_Coding_Region,Number_of_Matched_ORFs_Extending_Start_Region,Percentage_of_Matched_ORFs_Extending_Start_Region,Number_of_Matched_ORFs_Extending_Stop_Region,Percentage_of_Matched_ORFs_Extending_Stop_Region,Number_of_All_ORFs_on_Positive_Strand,Percentage_of_All_ORFs_on_Positive_Strand,Number_of_All_ORFs_on_Negative_Strand,Percentage_of_All_ORFs_on_Negative_Strand,Median_Start_Difference_of_Matched_ORFs,Median_Stop_Difference_of_Matched_ORFs,ATG_Start_Percentage,GTG_Start_Percentage,TTG_Start_Percentage,ATT_Start_Percentage,CTG_Start_Percentage,Other_Start_Codon_Percentage,TAG_Stop_Percentage,TAA_Stop_Percentage,TGA_Stop_Percentage,Other_Stop_Codon_Percentage,True_Positive,False_Positive,False_Negative,Precision,Recall,False_Discovery_Rate,Nucleotide_True_Positive,Nucleotide_False_Positive,Nucleotide_True_Negative,Nucleotide_False_Negative,Nucleotide_Precision,Nucleotide_Recall,Nucleotide_False_Discovery_Rate,ORF_Nucleotide_Coverage_of_Genome,Matched_ORF_Nucleotide_Coverage_of_Genome
4263,5.21,3971,93.15,3973,98.05,824.0,-2.83,89,102.27,7103,0.38,52.06,-0.18,52.21,0.11,1592,74.37,143,0.00,878,-3.83,112,3.00,468,12.77,346,-16.63,3737,94.11,3737,94.11,3973,100.05,2,0,0.00,97,2.44,0,0.00,2074,0.49,2189,0.51,-3.0,N/A,90.62,7.65,1.71,0.00,0.00,0.02,7.79,63.31,28.90,0.00,0.98,0.07,0.02,0.93,0.98,0.07,1.00,0.23,0.77,0.00,0.96,1.00,0.04,87.56,83.96
Reference_CDS_Gene_Coverage_of_Genome
84.35
Predicted_CDS_Coverage_of_Genome
87.56
Matched_Predicted_CDS_Coverage_of_Genome
83.96
Start_Position_Difference:
57,6,-126,-111,-75,-3,-39,-33,3,108,6,21,-138,-9,27,123,78,414,-36,-24,-126,-12,-126,-75,-48,3,-39,-27,-6,-66,9,66,84,-3,33,3,-78,-54,-39,-33,6,3,-36,-126,180,-99,123,-78,-30,72,-36,39,69,-15,60,-3,-33,-135,-3,-6,-81,-48,21,108,9,-15,126,-45,3,-57,-30,-60,33,30,-27,3,39,30,63,-3,48,6,-111,30,-60,15,-27,66,-21,-39,-45,-15,-60,12,36,-123,33,45,-24,36,-21,-57,636,-12,-9,36,27,-12,-3,-6,-48,-9,-108,21,-3,-3,15,-12,66,-3,54,-84,54,57,57,-63,111,216,-57,-27,-33,3,-3,-36,-54,21,-36,3,-81,-45,-30,-126,-39,-36,-6,-18,3,12,42,-3,-15,21,45,18,3,-39,825,-45,51,3,-75,45,-27,105,3,72,30,189,60,-39,3,3,-141,6,-54,42,-72,36,-15,27,-111,120,-90,6,54,-108,9,15,-36,-3,12,-30,126,33,-36,-72,12,-72,-39,-36,12,-126,15,36,-3,48,72,-33,33,21,-6,48,-6,-24,-636,-90,6,-105,-57,-24,60,-126,54,-45,36,-468,9,3,12,45,12,15,36,3,21,30
Stop_Position_Difference:

Alternative_Starts_Predicted:
CTT
Alternative_Stops_Predicted:

Undetected_Gene_Metrics:
ATG_Start ,GTG_Start ,TTG_Start ,ATT_Start ,CTG_Start ,Alternative_Start_Codon ,TGA_Stop ,TAA_Stop ,TAG_Stop ,Alternative_Stop_Codon ,Median_Length ,ORFs_on_Positive_Strand ,ORFs_on_Negative_Strand
86.08,10.13,3.80,0.00,0.00,0.00,34.18,56.96,7.59,1.27,128.00,39,40
Perfect_Match_Genes:
>E-coli_337_2799_+ 
ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGATGGCAGGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAGCTTTCCTACTTCGGCGCTAAAGTTCTTCACCCCCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCTGATTAAAAATACCGGAAATCCTCAAGCACCAGGTACGCTCATTGGTGCCAGCCGTGATGAAGACGAATTACCGGTCAAGGGCATTTCCAATCTGAATAACATGGCAATGTTCAGCGTTTCTGGTCCGGGGATGAAAGGGATGGTCGGCATGGCGGCGCGCGTCTTTGCAGCGATGTCACGCGCCCGTATTTCCGTGGTGCTGATTACGCAATCATCTTCCGAATACAGCATCAGTTTCTGCGTTCCACAAAGCGACTGTGTGCGAGCTGAACGGGCAATGCAGGAAGAGTTCTACCTGGAACTGAAAGAAGGCTTACTGGAGCCGCTGGCAGTGACGGAACGGCTGGCCATTATCTCGGTGGTAGGTGATGGTATGCGCACCTTGCGTGGGATCTCGGCGAAATTCTTTGCCGCACTGGCCCGCGCCAATATCAACATTGTCGCCATTGCTCAGGGATCTTCTGAACGCTCAATCTCTGTCGTGGTAAATAACGATGATGCGACCACTGGCGTGCGCGTTACTCATCAGATGCTGTTCAATACCGATCAGGTTATCGAAGTGTTTGTGATTGGCGTCGGTGGCGTTGGCGGTGCGCTGCTGGAGCAACTGAAGCGTCAGCAAAGCTGGCTGAAGAATAAACATATCGACTTACGTGTCTGCGGTGTTGCCAACTCGAAGGCTCTGCTCACCAATGTACATGGCCTTAATCTGGAAAACTGGCAGGAAGAACTGGCGCAAGCCAAAGAGCCGTTTAATCTCGGGCGCTTAATTCGCCTCGTGAAAGAATATCATCTGCTGAACCCGGTCATTGTTGACTGCACTTCCAGCCAGGCAGTGGCGGATCAATATGCCGACTTCCTGCGCGAAGGTTTCCACGTTGTCACGCCGAACAAAAAGGCCAACACCTCGTCGATGGATTACTACCATCAGTTGCGTTATGCGGCGGAAAAATCGCGGCGTAAATTCCTCTATGACACCAACGTTGGGGCTGGATTACCGGTTATTGAGAACCTGCAAAATCTGCTCAATGCAGGTGATGAATTGATGAAGTTCTCCGGCATTCTTTCTGGTTCGCTTTCTTATATCTTCGGCAAGTTAGACGAAGGCATGAGTTTCTCCGAGGCGACCACGCTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTTTCTGGTATGGATGTGGCGCGTAAACTATTGATTCTCGCTCGTGAAACGGGACGTGAACTGGAGCTGGCGGATATTGAAATTGAACCTGTGCTGCCCGCAGAGTTTAACGCCGAGGGTGATGTTGCCGCTTTTATGGCGAATCTGTCACAACTCGACGATCTCTTTGCCGCGCGCGTGGCGAAGGCCCGTGATGAAGGAAAAGTTTTGCGCTATGTTGGCAATATTGATGAAGATGGCGTCTGCCGCGTGAAGATTGCCGAAGTGGATGGTAATGATCCGCTGTTCAAAGTGAAAAATGGCGAAAACGCCCTGGCCTTCTATAGCCACTATTATCAGCCGCTGCCGTTGGTACTGCGCGGATATGGTGCGGGCAATGACGTTACAGCTGCCGGTGTCTTTGCTGATCTGCTACGTACCCTCTCATGGAAGTTAGGAGTCTGA 
Partial_Match_Genes:
Gene:16751_16903_-_ATG_TAA 
ATGAAGCAGCATAAGGCGATGATTGTCGCCCTGATCGTCATCTGTATCACCGCCGTAGTGGCGGCGCTGGTAACGAGAAAAGACCTCTGTGAGGTTCACATCCGAACTGGCCAGACGGAGGTTGCTGTTTTCACGGCTTACGAATCCGAGTAA 
Predicted_CDS:16751_16960_-_ATG_TAA 
ATGCTGAACACATGTAGAGTGCCTCTTACTGACCGTAAGGTCAAGGAGAAGAGAGCAATGAAGCAGCATAAGGCGATGATTGTCGCCCTGATCGTCATCTGTATCACCGCCGTAGTGGCGGCGCTGGTAACGAGAAAAGACCTCTGTGAGGTTCACATCCGAACTGGCCAGACGGAGGTTGCTGTTTTCACGGCTTACGAATCCGAGTAA 
Missed_Genes:
>E-coli_190_255_+ 
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA 
Predicted_CDSs_Without_Corresponding_Gene_In_Reference_Metrics:
ATG_Start ,GTG_Start ,TTG_Start ,ATT_Start ,CTG_Start ,Alternative_Start_Codon ,TGA_Stop ,TAA_Stop ,TAG_Stop ,Alternative_Stop_Codon ,Median_Length ,ORFs_on_Positive_Strand ,ORFs_on_Negative_Strand
82.88,12.33,4.45,0.00,0.00,0.34,36.99,42.81,20.21,0.00,356.00,152,140
Predicted_CDS_Without_Corresponding_Gene_in_Reference:
>Prodigal_3_98_+ 
CTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAA
Predicted_CDSs_Which_Detected_more_than_one_Gene:
Predicted_CDS:16751-16960_Genes:16751-16960|16751-16903
Predicted_CDS:2652100-2652990_Genes:2652100-2652990|2652100-2652165
```


## GFF Tools:

### GFF_Adder:

GFF_Adder allows for the addition of predicted CDSs to an existing reference annotation (GFF or another tool) which produces a new GFF containing the original
genes plus the new CDS from another prediction. Default filtering will remove additional CDSs that overlap existing genes by more than 50 nt.
The ```-gi``` option can be used to allow for different genomic elements to be accounted for, other than only CDSs in the reference annotation.

For Help: ```python python3 -m ORForise.GFF_Adder -h ```

```python
usage: GFF_Adder.py [-h] -dna GENOME_DNA [-rt REFERENCE_TOOL] -ref REFERENCE_ANNOTATION [-gi GENE_IDENT] -at ADDITIONAL_TOOL -add ADDITIONAL_ANNOTATION [-olap OVERLAP] -o OUTPUT_FILE

optional arguments:
  -h, --help            show this help message and exit
  -dna GENOME_DNA, --genome_DNA GENOME_DNA
                        Genome DNA file (.fa) which both annotations are based on
  -rt REFERENCE_TOOL, --reference_tool REFERENCE_TOOL
                        Which tool format to use as reference? - If not provided, will default to standard Ensembl GFF format, can be Prodigal or any of the other tools available
  -ref REFERENCE_ANNOTATION, --reference_annotation REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -gi GENE_IDENT, --gene_ident GENE_IDENT
                        Identifier used for extraction of "genic" regions from reference annotation "CDS,rRNA,tRNA": Default for is "CDS"
  -at ADDITIONAL_TOOL, --additional_tool ADDITIONAL_TOOL
                        Which format to use for additional annotation?
  -add ADDITIONAL_ANNOTATION, --additional_annotation ADDITIONAL_ANNOTATION
                        Which annotation file to add to reference annotation?
  -olap OVERLAP, --overlap OVERLAP
                        Maximum overlap between reference and additional genic regions (CDS,rRNA etc) - Default: 50 nt
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output filename

```

### GFF_Intersector:

GFF_Intersector enables the aggregation of different genome annotations and CDS predictions and creates a single GFF
representing the intersection of the two existing annotations.
GFF_Intersector also provides an option to allow the retention of genes that have a user defined difference (minimum % coverage and in-frame).
The ```-gi``` option can be used to allow for different genomic elements to be accounted for, other than only CDSs in the reference annotation.

For Help: ```python python3 -m ORForise.GFF_Intersector -h ``` 
```python
usage: GFF_Intersector.py [-h] -dna GENOME_DNA [-rt REFERENCE_TOOL] -ref REFERENCE_ANNOTATION [-gi GENE_IDENT] -at ADDITIONAL_TOOL -add ADDITIONAL_ANNOTATION [-cov COVERAGE] -o OUTPUT_FILE

optional arguments:
  -h, --help            show this help message and exit
  -dna GENOME_DNA, --genome_DNA GENOME_DNA
                        Genome DNA file (.fa) which both annotations are based on
  -rt REFERENCE_TOOL, --reference_tool REFERENCE_TOOL
                        Which tool format to use as reference? - If not provided, will default to standard Ensembl GFF format, can be Prodigal or any of the other tools available
  -ref REFERENCE_ANNOTATION, --reference_annotation REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -gi GENE_IDENT, --gene_ident GENE_IDENT
                        Identifier used for extraction of "genic" regions from reference annotation "CDS,rRNA,tRNA": Default for is "CDS"
  -at ADDITIONAL_TOOL, --additional_tool ADDITIONAL_TOOL
                        Which format to use for additional annotation?
  -add ADDITIONAL_ANNOTATION, --additional_annotation ADDITIONAL_ANNOTATION
                        Which annotation file to add to reference annotation?
  -cov COVERAGE, --coverage COVERAGE
                        Percentage coverage of reference annotation needed to confirm intersection - Default: 100 == exact match
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output filename

```

# Genomes Available:

The .fa and .gff files (from Ensembl Bacteria) below are available in the Genomes directory.

* *Escherichia coli K-12* - Strain ER3413 - Assembly ASM80076v1
* *Staphylococcus aureus* - Strain 502A - Assembly ASM59796v1
* *Bacillus subtilis* - Strain BEST7003 - Assembly ASM52304v1
* *Mycoplasma genitalium* - Strain G37 - Assembly ASM2732v1
* *Caulobacter crescentus* - Strain CB15 - Assembly ASM690v1
* *Pseudomonas fluorescens* - Strain UK4 - Assembly ASM73042v1

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


