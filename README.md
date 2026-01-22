# ORForise - Genome Annotation Analysis and Comparison Platform
## Published in Bioinformatics :   https://academic.oup.com/bioinformatics/article/38/5/1198/6454948

# Requirements and Installation:

### The ORForise platform is written in Python (3.6-3.*) and only requires the NumPy library (should be installed automatically by pip when installing ORForise) which is standard in most base installations of Python3.

## Intallation:

### ORForise is available via the pip Python package manager ```pip3 install ORForise``` and bioconda ```conda install -c bioconda ORForise```.

### Testing:
Precomputed testing and data which includes example input and output files for all tools presented below is available in the `~ORForise/Testing` directory of the GitHub repository. 
Example output files from ```Annotation-Compare```, ```Aggregate-Compare```, ```Convert-To-GFF``` and ```Annotation-Intersector``` are available.


## Genome Annotation Analysis:

### Use-cases: (Running if via pip)

For Help: ```Annotation-Compare -h ```

```python
ORForise v1.6.4: Annotatione-Compare Run Parameters.

Required Arguments:
  -dna GENOME_DNA       Genome DNA file (.fa) which both annotations are based on
  -ref REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -t TOOL               Which tool to analyse? 
  -tp TOOL_PREDICTION   Tool genome prediction file (.gff) - Different Tool Parameters are compared individually via separate files

Optional Arguments:
  -gene_ident GENE_IDENT
                        What features to consider as genes? - Default: CDS - Provide comma separated list of features to consider as genes (e.g. CDS,exon)
  -rt REFERENCE_TOOL    What type of Annotation to compare to? -- Leave blank for Ensembl reference- Provide tool name to compare output from two tools

Output:
  -o OUTDIR             Define directory where detailed output should be places
  -n OUTNAME            Define output filename(s) prefix - If not provided, filename of reference annotation file will be used- <outname>_<contig_id>_ORF_Comparison.csv

Misc:
  -v {True,False}       Default - False: Print out runtime status
```
### Compare a *de novo* genome annotation to an Ensembl annotation:

Genome annotation is a difficult process, even for Prokaryotes. ORForise allows for the direct and systematic analysis of
*de novo* gene prediction from a wide selection of tools to a reference Genome Annotation, such as those provided by
Ensembl Bacteria.

#### Example: Installation through pip will allow user to call the programs directly from the ORForise package (Prodigal and Pyrodigal provide annotations in the same format).
```python
Annotation-Compare -dna ~/Test_Data/Genomes/E-coli/Escherichia_coli.fasta -ref ~/Test_Data/Genomes/E-coli/Escherichia_coli.gff -t Prodigal -tp ~/Test_Data/Genomes/E-coli/Prodigal_Escherichia_coli.gff
```
### Example Output: - See ```~/Test_Data/Genomes/E-coli/annotation_compare```
```commandline
Genome Used: Escherichia_coli.fasta
Reference Used: Escherichia_coli.gff
Tool Compared: Prodigal
Total Number of Reference Genes: 5222
Number of Contigs: 4
Contig	Genes	ORFs	Perfect_Matches	Partial_Matches	Missed_Genes	Unmatched_ORFs	Multi_Matched_ORFs
ERS715463SCcontig000003	4068	4070	4065	1	2	4	0
ERS715463SCcontig000002	1033	1035	1033	0	0	2	0
ERS715463SCcontig000001	75	77	75	0	0	2	0
ERS715463SCcontig000004	46	47	45	1	0	1	0

Overall Summary:
Number of Genes: 5222
Number of ORFs: 5229
Perfect Matches: 5218 [5222] - 99.92%
Partial Matches: 2 [5222] - 0.04%
Missed Genes: 2 [5222] - 0.04%
Unmatched ORFs: 9 [5222] - 0.17%
Multi-matched ORFs: 0 [5222] - 0.00%

```


## Compare different novel annotations with each other on a single Genome:

If a reference Genome Annotation is not available or a direct comparison between two or more tools is wanted,
ORForise can be used as the example below.

## Aggregate CDS Prediction Analysis:

### Use-cases: (Running if via pip)

For Help: ```Aggregate-Compare -h ```

```python
ORForise v1.6.4: Aggregate-Compare Run Parameters.

Required Arguments:
  -dna GENOME_DNA       Genome DNA file (.fa) which both annotations are based on
  -t TOOLS              Which tools to analyse?
  -tp TOOL_PREDICTIONS  Tool genome prediction file (.gff) - Providefile locations for each tool comma separated
  -ref REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?

Optional Arguments:
  -gene_ident GENE_IDENT
                        What features to consider as genes? - Default: CDS - Provide comma separated list of features to consider as genes (e.g. CDS,exon)
  -rt REFERENCE_TOOL    What type of Annotation to compare to? -- Leave blank for Ensembl reference- Provide tool name to compare output from two tools

Output:
  -o OUTDIR             Define directory where detailed output should be places - If not provided, summary will be printed to std-out
  -n OUTNAME            Define output file name - Mandatory is -o is provided: <outname>_<contig_id>_ORF_Comparison.csv

Misc:
  -v {True,False}       Default - False: Print out runtime status

```

#### Example: 
```python
Aggregate-Compare -ref ~/Test_Data/Genomes/E-coli/Escherichia_coli.gff -dna ~/Test_Data/Genomes/E-coli/Escherichia_coli.fasta -t Prodigal,GeneMarkS2 -tp ~/Test_Data/Genomes/E-coli/Prodigal_Escherichia_coli.gff,~/Test_Data/Genomes/E-coli/GeneMarkS2_E-coli.gff
```
This will compare and agregate the predictions of Prodigal and GeneMarkS2 against the E-coli reference annotation provided by Ensembl Bacteria.

### Annotation Comparison Output - The output format is the same for Annotation_Compare and Aggregate_Compare:  See ```~/Test_Data/Genomes/E-coli/aggregate_compare```
```bash
Genome Used: Escherichia_coli.fasta
Reference Used: Escherichia_coli.gff
Tool Compared: Prodigal,GeneMarkS2
Total Number of Reference Genes: 5222
Number of Contigs: 4
Contig	Genes	ORFs	Perfect_Matches	Partial_Matches	Missed_Genes	Unmatched_ORFs	Multi_Matched_ORFs
ERS715463SCcontig000003	4068	4500	4065	1	2	434	0
ERS715463SCcontig000002	1033	1148	1033	0	0	115	0
ERS715463SCcontig000001	75	92	75	0	0	17	0
ERS715463SCcontig000004	46	64	45	1	0	18	0

Overall Summary:
Number of Genes: 5222
Number of ORFs: 5804
Perfect Matches: 5218 [5222] - 99.92%
Partial Matches: 2 [5222] - 0.04%
Missed Genes: 2 [5222] - 0.04%
Unmatched ORFs: 584 [5222] - 11.18%
Multi-matched ORFs: 0 [5222] - 0.00%

  Prodigal: Perfect=5218, Partial=2, Unmatched=9, Multi-matched=0

  GeneMarkS2: Perfect=4609, Partial=2, Unmatched=579, Multi-matched=0
```

Shown so far have been the summary outputs of the comparison tools.
Since v 1.5.0, detailed CSV outputs are also provided for each contig analysed - See ```~/Test_Data/Genomes/E-coli/annotation_compare``` for example outputs.
```commandline
Representative_Metrics:
Percentage_of_Genes_Detected,Percentage_of_ORFs_that_Detected_a_Gene,Percent_Difference_of_All_ORFs,Median_Length_Difference,Percentage_of_Perfect_Matches,Median_Start_Difference_of_Matched_ORFs,Median_Stop_Difference_of_Matched_ORFs,Percentage_Difference_of_Matched_Overlapping_CDSs,Percent_Difference_of_Short-Matched-ORFs,Precision,Recall,False_Discovery_Rate
100.00,97.87,2.17,1.20,97.83,6.0,N/A,0.00,0.00,0.98,1.00,0.02
Prediction_Metrics:
Number_of_ORFs,Percent_Difference_of_All_ORFs,Number_of_ORFs_that_Detected_a_Gene,Percentage_of_ORFs_that_Detected_a_Gene,Number_of_Genes_Detected,Percentage_of_Genes_Detected,Median_Length_of_All_ORFs,Median_Length_Difference,Minimum_Length_of_All_ORFs,Minimum_Length_Difference,Maximum_Length_of_All_ORFs,Maximum_Length_Difference,Median_GC_content_of_All_ORFs,Percent_Difference_of_All_ORFs_Median_GC,Median_GC_content_of_Matched_ORFs,Percent_Difference_of_Matched_ORF_GC,Number_of_ORFs_which_Overlap_Another_ORF,Percent_Difference_of_Overlapping_ORFs,Maximum_ORF_Overlap,Median_ORF_Overlap,Number_of_Matched_ORFs_Overlapping_Another_ORF,Percentage_Difference_of_Matched_Overlapping_CDSs,Maximum_Matched_ORF_Overlap,Median_Matched_ORF_Overlap,Number_of_Short-ORFs,Percent_Difference_of_Short-ORFs,Number_of_Short-Matched-ORFs,Percent_Difference_of_Short-Matched-ORFs,Number_of_Perfect_Matches,Percentage_of_Perfect_Matches,Number_of_Perfect_Starts,Percentage_of_Perfect_Starts,Number_of_Perfect_Stops,Percentage_of_Perfect_Stops,Number_of_Out_of_Frame_ORFs,Number_of_Matched_ORFs_Extending_a_Coding_Region,Percentage_of_Matched_ORFs_Extending_a_Coding_Region,Number_of_Matched_ORFs_Extending_Start_Region,Percentage_of_Matched_ORFs_Extending_Start_Region,Number_of_Matched_ORFs_Extending_Stop_Region,Percentage_of_Matched_ORFs_Extending_Stop_Region,Number_of_All_ORFs_on_Positive_Strand,Percentage_of_All_ORFs_on_Positive_Strand,Number_of_All_ORFs_on_Negative_Strand,Percentage_of_All_ORFs_on_Negative_Strand,Median_Start_Difference_of_Matched_ORFs,Median_Stop_Difference_of_Matched_ORFs,ATG_Start_Percentage,GTG_Start_Percentage,TTG_Start_Percentage,ATT_Start_Percentage,CTG_Start_Percentage,Other_Start_Codon_Percentage,TAG_Stop_Percentage,TAA_Stop_Percentage,TGA_Stop_Percentage,Other_Stop_Codon_Percentage,True_Positive,False_Positive,False_Negative,Precision,Recall,False_Discovery_Rate,Nucleotide_True_Positive,Nucleotide_False_Positive,Nucleotide_True_Negative,Nucleotide_False_Negative,Nucleotide_Precision,Nucleotide_Recall,Nucleotide_False_Discovery_Rate,ORF_Nucleotide_Coverage_of_Genome,Matched_ORF_Nucleotide_Coverage_of_Genome
47,2.17,46,97.87,46,100.00,378.0,1.20,138,0.00,1551,0.00,54.79,-0.05,54.81,0.00,36,63.64,61,0.00,22,0.00,61,3.00,13,0.00,13,0.00,45,97.83,45,97.83,46,100.00,0,0,0.00,1,2.17,0,0.00,11,0.23,36,0.77,6.0,N/A,82.98,12.77,2.13,0.00,0.00,2.13,6.38,23.40,68.09,2.13,1.00,0.02,0.00,0.98,1.00,0.02,1.00,0.16,0.84,0.00,0.96,1.00,0.04,81.98,78.63
Reference_CDS_Gene_Coverage_of_Genome
78.61
Predicted_CDS_Coverage_of_Genome
81.98
Matched_Predicted_CDS_Coverage_of_Genome
78.63

```


## GFF/Annotation Manipulation Tools: ORForise also provides tools to manipulate and combine existing annotations in GFF format or other tool-specific formats.
### GFF-Adder:
GFF-Adder combines two existing annotations (GFF or other tool formats).
For Help: ```GFF-Adder -h ```

```python
ORForise v1.6.4: GFF-Adder Run Parameters.

Required Arguments:
  -dna GENOME_DNA       Genome DNA file (.fa) which both annotations are based on
  -ref REFERENCE_ANNOTATION
                        Which reference annotation file to use as reference?
  -at ADDITIONAL_TOOL   Which format to use for additional annotation? - Can provide multiple annotations (Tool1,Tool2)
  -add ADDITIONAL_ANNOTATION
                        Which annotation file to add to reference annotation? - Can provide multiple annotations (1.GFF,2.GFF)
  -o OUTPUT_FILE        Output filename

Optional Arguments:
  -rt REFERENCE_TOOL    Which tool format to use as reference? - If not provided, will default to the standard GFF format and will only look for "CDS" features
  --gene_ident GENE_IDENT
                        Identifier used for identifying genomic features in reference annotation "CDS,rRNA,tRNA"
  -mc                   Default - False: Mark reference annotations which where present in the additional tool annotation
  -c                    Default - False: Do not mark 9th column with "Original/Matched/Additional tag"
  --meta                Default - False: Output metadata file
  --olap OVERLAP        Maximum overlap between reference and additional genic regions (CDS,rRNA etc) - Default: 50 nt

Misc:
  -v {True,False}       Default - False: Print out runtime status


```

#### Example: Running GFF-Adder to combine the additional CDS predictions made by Prodial to the canonical annotations from Ensembl.
``` GFF-Adder -dna ~/Test_Data/Genomes/E-coli/Escherichia_coli.fasta -ref ~/Test_Data/Genomes/E-coli/Escherichia_coli.gff  -at Prodigal -add ~/Test_Data/Genomes/E-coli/Prodigal_Escherichia_coli.gff -o ~/Test_Data/Genomes/E-coli/Ensembl_AND_Prodigal_Escherichia_coli.gff ```
#### Example Output: ~/ORForise/Testing/Myco_Ensembl_GFF_Adder_Prodigal.gff
```
##gff-version	3
#	GFF-Adder
#	Run Date:2026-01-11
##Genome DNA File:../../Test_Data/Genomes/E-coli/Escherichia_coli.fasta
##Original File: ../../Test_Data/Genomes/E-coli/Escherichia_coli.gff
##Additional File: ../../Test_Data/Genomes/E-coli/Prodigal_Escherichia_coli.gff
ERS715463SCcontig000003	Prodigal	CDS	2	388	.	+	.	ID=Additional_Annotations;Prodigal
ERS715463SCcontig000003	MGnify	CDS	83	388	.	+	.	ID=Original_Annotation;ID=ENSB_0kRwXBh8bjHtVl3;Parent=transcript:ENSB:0kRwXBh8bjHtVl3;protein_id=ENSB:0kRwXBh8bjHtVl3
ERS715463SCcontig000003	MGnify	CDS	453	542	.	+	.	ID=Original_Annotation;ID=ENSB_W8Go0tx9y9dAtng;Parent=transcript:ENSB:W8Go0tx9y9dAtng;protein_id=ENSB:W8Go0tx9y9dAtng;Matched_Annotations=Prodigal
```

### Annotation-Intersector:

Annotation-Intersector combines and contracts two existing annotations (GFF or other tool formats)

For Help: ```Annotation-Intersector -h ``` 
```python
Thank you for using ORForise
Please report any issues to: https://github.com/NickJD/ORForise/issues
#####
usage: Annotation_Intersector.py [-h] -ref REFERENCE_ANNOTATION -at
                                 ADDITIONAL_TOOL -add ADDITIONAL_ANNOTATION -o
                                 OUTPUT_FILE [-dna GENOME_DNA]
                                 [-rt REFERENCE_TOOL] [-gi GENE_IDENT]
                                 [-cov COVERAGE] [--report-discordance]
                                 [--report-discordance-file REPORT_DISCORDANCE_FILE]

ORForise v1.6.4: Annotation-Intersector Run Parameters

options:
  -h, --help            show this help message and exit

Required Arguments:
  -ref REFERENCE_ANNOTATION
                        Reference annotation GFF file
  -at ADDITIONAL_TOOL   Tool name/format for additional annotation (module
                        under Tools/)
  -add ADDITIONAL_ANNOTATION
                        Additional annotation file to compare
  -o OUTPUT_FILE        Output GFF filename for kept genes

Optional Arguments:
  -dna GENOME_DNA       Genome DNA file (.fa) which both annotations are based
                        on
  -rt REFERENCE_TOOL    Reference tool parser name (if not provided, GFF is
                        expected)
  -gi GENE_IDENT        Comma-separated feature types to consider from
                        reference (default: CDS)
  -cov COVERAGE, --coverage COVERAGE
                        Percentage coverage threshold for intersection
                        (default 100)
  --report-discordance  If set, produce discordance reports (three GFFs)
  --report-discordance-file REPORT_DISCORDANCE_FILE
                        Optional base path for discordance reports

```

#### Example: Running Annotation-Intersector to combine and contract annotations from multiple tools or reference files.
``` Annotation-Intersector  -ref .../ORForise/Tools/EasyGene/EasyGene_E-coli_E-coli.gff -rt EasyGene -at Prodigal -add .../ORForise/Tools/Prodigal/Prodigal_E-coli.gff -o .../Test_Data/Annotation-Intersector/Annotation-Intersect.gff --report-discordance ```

#### Example Output:
##### .../Test_Data/Annotation-Intersector/Annotation-Intersect.gff
```
##gff-version	3
#	Annotation-Intersector
#	Run Date:2026-01-09
##Original File: .../ORForise/Tools/EasyGene/EasyGene_E-coli_E-coli.gff
##Intersecting File: .../ORForise/Tools/Prodigal/Prodigal_E-coli.gff
Chromosome	EasyGene	CDS	337	2799	.	+	.	ID=Original_Annotation=EasyGene;Additional_Annotation=Prodigal;Coverage=100.0
Chromosome	EasyGene	CDS	3734	5020	.	+	.	ID=Original_Annotation=EasyGene;Additional_Annotation=Prodigal;Coverage=100.0
Chromosome	EasyGene	CDS	5683	6459	.	-	.	ID=Original_Annotation=EasyGene;Additional_Annotation=Prodigal;Coverage=100.0
Chromosome	EasyGene	CDS	6529	7959	.	-	.	ID=Original_Annotation=EasyGene;Additional_Annotation=Prodigal;Coverage=100.0
Chromosome	EasyGene	CDS	8238	9191	.	+	.	ID=Original_Annotation=EasyGene;Additional_Annotation=Prodigal;Coverage=100.0
```
#### .../Test_Data/Annotation-Intersector/Annotation-Intersect.only_in_reference.gff
```
##gff-version	3
#	Annotation-Intersector discordance report
#	Run Date:2026-01-09
##Original File: EasyGene_E-coli_E-coli
Chromosome	EasyGene	CDS	408401	408484	.	.	.	Status=only_in_ref;Coverage=0.00;Ref_info=EasyGene
Chromosome	EasyGene	CDS	1272584	1272886	.	.	.	Status=only_in_ref;Coverage=0.00;Ref_info=EasyGene
Chromosome	EasyGene	CDS	2574901	2574960	.	.	.	Status=only_in_ref;Coverage=0.00;Ref_info=EasyGene
Chromosome	EasyGene	CDS	2710019	2710081	.	.	.	Status=only_in_ref;Coverage=0.00;Ref_info=EasyGene
```
#### .../Test_Data/Annotation-Intersector/Annotation-Intersect.mismatches.gff
```
##gff-version	3
#	Annotation-Intersector discordance report
#	Run Date:2026-01-09
##Original File: EasyGene_E-coli_E-coli
Chromosome	EasyGene	CDS	18715	19620	.	.	.	Status=found_in_additional_but_below_coverage;Coverage=99.34;Ref_info=EasyGene;Add_info=Prodigal
Chromosome	EasyGene	CDS	19811	20314	.	.	.	Status=found_in_additional_but_below_coverage;Coverage=75.00;Ref_info=EasyGene;Add_info=Prodigal
Chromosome	EasyGene	CDS	29624	30799	.	.	.	Status=found_in_additional_but_below_coverage;Coverage=97.70;Ref_info=EasyGene;Add_info=Prodigal
Chromosome	EasyGene	CDS	70378	71265	.	.	.	Status=found_in_additional_but_below_coverage;Coverage=98.99;Ref_info=EasyGene;Add_info=Prodigal

```

#### Convert-To-GFF: Converts tool-specific output files to standard GFF3 format for use in ORForise analyses.
For Help: ```Convert_To_GFF.py -h ```
```
Thank you for using ORForise
Please report any issues to: https://github.com/NickJD/ORForise/issues
#####
usage: Convert_To_GFF.py [-h] [-dna GENOME_DNA] -i INPUT_ANNOTATION -fmt FORMAT -o OUTPUT_DIR [-gi GENE_IDENT] [--verbose]

ORForise v1.6.4: Convert-To-GFF Run Parameters

Required Arguments:
  -dna GENOME_DNA      Genome DNA file (.fa)
  -i INPUT_ANNOTATION  Input annotation file (tabular)
  -fmt FORMAT          Input format: blast, abricate, genemark
  -o OUTPUT_DIR        Output directory

Optional Arguments:
  -gi GENE_IDENT       Gene identifier types to extract (unused)
  --verbose            Verbose logging with logfile
```

# Genomes Available:

The .fa and .gff files (from Ensembl Bacteria Release 46) below are available in the Genomes directory.

* *Bacillus subtilis* - Strain BEST7003 - Assembly ASM52304v1
* *Caulobacter crescentus* - Strain CB15 - Assembly ASM690v1
* *Escherichia coli K-12* - Strain ER3413 - Assembly ASM80076v1
* *Mycoplasma genitalium* - Strain G37 - Assembly ASM2732v1
* *Pseudomonas fluorescens* - Strain UK4 - Assembly ASM73042v1
* *Staphylococcus aureus* - Strain 502A - Assembly ASM59796v1



# Prediction Tool Formats Currently Available:
ORForise currently supports the comparison of multiple gene prediction tools via their output in GFF3 format. \
This can be used to compare different annotations with eachother or additional tools which use the GFF3 format.

## Tool Specific Formats:
Run ```List-Tools``` to see the available tools. \
ORForise only needs the tool name and the annotation file produced from any compatible tool to undertake the analysis.

**If the tool uses another non-standard format, a request can be made to add it as an option via GitHub.**

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

**GeneMarkHA - Version 3.25** - http://exon.gatech.edu/GeneMark/heuristic_gmhmmp.cgi  
GFF was chosen as output type.

**GeneMarkS - Version 4.25** - http://exon.gatech.edu/GeneMark/genemarks.cgi  
GFF was chosen as output type.

**GeneMarkS2 - Version '2020'** - http://exon.gatech.edu/GeneMark/genemarks2.cgi   
GFF3 was chosen as output type.

**GLIMMER3 - Version 3.02** - http://ccb.jhu.edu/software/glimmer/index.shtml  
Default parameters from manual were used.

**MetaGene - Version 2.24.0** - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1636498/  
Default options were used.

**MetaGeneAnnotator - Version '2008/8/19'** - http://metagene.nig.ac.jp/  
Defaults options were used.

**MetaGeneMark - Version '2020'** - http://exon.gatech.edu/meta_gmhmmp.cgi  
GFF was chosen as output type.

**Prodigal (Includes Pyrodigal) - Version 2.6.3** - https://github.com/hyattpd/Prodigal  
GFF was chosen as output type.

**TransDecoder - Version 5.5.0** - https://github.com/TransDecoder/TransDecoder/wiki  
Defaults options were used.

**Balrog - Version 2021`** - https://github.com/salzberg-lab/Balrog  
Defaults options were used.


