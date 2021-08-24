# ORForise - Prokaryote Genome Annotation Comparison and Analysis Platform

Platform for analysing Prokaryote CoDing Sequence (CDS) Gene Predictors. \
Novel genome annotations can be compared to a provided reference annotation from Ensembl (or any given GFF annotation) 
and predictions from other tools.

## Requirements and Installation:

The ORForise platform is written in Python3.8 and only requires the NumPy library which is standard in most base
installations of Python3. \
Usually, ```pip3 install numpy``` is adequate to install NumPy.

### Intallation:

The ORForise platform is available via github ```git clone https://github.com/NickJD/ORForise``` and the pip Python package manager ```pip install ORForise```. \
For both methods of 'installation', cloning the repository or installing ORForise via pip, the same input files are required when running the code. 

To run, you need:
* Input Genome FASTA and corresponding GFF file or CDS predictions with the annotated genes for the genome you want to use as reference.
* A prediction output from one of the compatible tools for the same genome.

Each tool requires its own directory and prediction output must be named as follows '~toolName/toolName_Species.*'

### How to add your own Genome:

Corresponding FASTA and GFF files must be provided for the genome the analysis is to be performed on.

### How to add your own tool:

If the new tool reports its predictions in GFF you can present ORForise with "GFF" for either the reference ```-rt``` or prediction ```-t``` option.
If the tool uses another non-standard format, a request can be made to add it as an option.

## CDS Prediction Analysis:

### Use-cases:

For Help: python3 -m ORForise.Annotation_Compare -h

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
a novel ORF prediction from a wide selection of tools to a reference Genome Annotation, such as those provided by
Ensembl Bacteria.

#### Example: Installation through pip will allow user to call the packages directly from the ORForise module.
```python
 ORForise.Annotation_Compare -dna Genomes/Myco.fa -ref Genomes/Myco.gff -t Prodigal -tp Tools/Prodigal/Prodigal_Myco.gff
```
### Compare different novel annotations with each other on a single Genome:

If a reference Genome Annotation is not available or a direct comparison between two or more tools is wanted,
ORForise can be used as the example below.

#### Example: 
```python
 ORForise.Annotation_Compare -dna Genomes/Myco.fa -rt  GeneMark_S_2 -ref Tools/GeneMark_S_2/GeneMark_S_2_Myco.gff -t Prodigal -tp Tools/Prodigal/Prodigal_Myco.gff
```
This will compare the novel Prodigal predictions against the predictions made by GeneMarkS-2

## GFF Tools:

### GFF_Adder:

GFF_Adder allows for the addition of predicted CDSs to an existing reference annotation (GFF or another tool) which produces a new GFF containing the original
genes plus the new CDS from another prediction. Default filtering will remove additional CDSs that overlap existing genes by more than 50 nt.
The ```-gi``` option can be used to allow for different genomic elements to be accounted for, other than only CDSs.

For Help: python3 -m ORForise.GFF_Adder -h

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
The ```-gi``` option can be used to allow for different genomic elements to be accounted for, other than only CDSs.

For Help: python3 -m ORForise.GFF_Intersector -h

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

The .fa and .gff files (from Ensembl Bacteria) for the Model Organisms below are available in the Genomes directory

* Escherichia coli K-12 - Strain ER3413 - Assembly ASM80076v1
* Staphylococcus aureus - Strain 502A - Assembly ASM59796v1
* Bacillus subtilis - Strain BEST7003 - Assembly ASM52304v1
* Mycoplasma genitalium - Strain G37 - Assembly ASM2732v1
* Caulobacter crescentus - Strain CB15 - Assembly ASM690v1
* Pseudomonas fluorescens - Strain UK4 - Assembly ASM73042v1

# Prediction Tools Available:

There are two Groups of tools - Those which do require a pre-built model and those which do not. \
They are listed with the non-default options used and their predictions for each of the 6 model organisms are available in their respective
directories.
ORForise only needs the tool name and the annotation file produced from any available model to undertake the analysis.

## GFF Standard Format:

The GFF Tool directory allows for the analysis of user-provided annotations in the standard GFF3 format.
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


