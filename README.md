# ORForise - Prokaryote Genome Annotation Analysis Platform 
Platform for analysing Prokaryote Protein Coding Gene Predictors.
Novel genome annotations can be compared to Gold Standard annotations from Ensembl (or any given GFF annotation) and predictions from other tools.


## Requirements and Installation:
ORForise is written in Python3 and only requires the Numpy library which is standard in most base installations of Python3. 
### Platform Preparation
First step is to clone the Repository `git clone https://github.com/NickJD/ORForise'

To run, you need:   
* Input Genome FASTA and corresponding GFF file with the annotated genes for the genome you want to compare to.  
* A prediction output from one of the compatible tools for the same genome (New tools can be added).   


Each tool requires it's own directory and prediction output must be named as follows '~toolName/toolName_Species.*'



## How to add your own Genome:
Each prediction tool has a directory and data handling script which converts the different tool output into a dictionary that
ORForise understands.  
New genomes and tools can be added to the analysis but they must follow the same directory and naming structure as originals.
To add a new genome, simply place the .fa and .gff files in the Genomes directory - Make sure the filename before the '.' is the same for both files.


## How to add your own tool:
If the new tool reports its predictions in GFF you can put the new tools output in the GFF directory. If another format is used, you must make another directory for that tool and reproduce the '~GFF/GFF.py' script so that ORForise knows how to extracting the ORF predictions.


# Use-cases:
For Help: python3 ORForise.py -h   
```python
usage: ORForise.py [-h] -t TOOL [-p PARAMETERS] -g GENOME_TO_COMPARE

arguments:
  -h, --help            show this help message and exit
  -t TOOL, --tool TOOL  Which tool to analyse?
  -p PARAMETERS, --parameters PARAMETERS
                        Optional parameters for prediction tool.
  -g GENOME_TO_COMPARE, --genome_to_compare GENOME_TO_COMPARE
                        Which genome to analyse? Genome files have same prefix - .fa and .gff appended
```
Output = "'tool'_'genome_to_compare'.csv"   

## Compare a novel genome annotation to an Ensembl Gold Standard: 
Genome annotation is a difficult process, 'even for Prokaryotes'. ORForise allows the direct and systematic analysis of a novel ORF prediction from a wide selection of tools to a Gold Standard Genome Annotation, such as those provided by Ensembl Bacteria.
### Example: python3 ORForise.py -t Prodigal -g E-coli 
### Example with Model Parameter: python3 ORForise.py -t GeneMark -p Staph -g Staph  

## Compare different novel annotations with eachother on a single Genome:
If a Gold Standard Genome Annotation is not available or a direct comparison between two or more tools is wanted, ORForise can be used as the example below.
### Example: python3 ORForise.py -t Prodigal -g E-coli_FragGeneScan
This will compare the novel Prodigal predictions against the predictions made by FragGeneScan (Note: For this example, the ~Genomes directory has both E-coli_FragGeneScan.fa and ~.gff files to allow the comparison)


# Genomes Compared:
The .fa and .gff files for the Model Organisms below are available in the ~Genomes directory
* Escherichia coli K-12 - Strain ER3413 - Assembly ASM80076v1    
* Staphylococcus aureus - Strain 502A - Assembly ASM59796v1  
* Bacillus subtilis - Strain BEST7003 - Assembly ASM52304v1  
* Mycoplasma genitalium - Strain G37 - Assembly ASM2732v1  
* Caulobacter crescentus - Strain CB15 - Assembly ASM690v1  
* Pseudomonas fluorescens - Strain UK4 - Assembly ASM73042v1  



# Prediction Tools Available:   
There are two Groups of tools - Those which do require a pre-built model and those which do not. 
They are listed with the non-default options used and their predictions for each of the 6 model organisms are available in their respective directories:
## Model Based Tools:

**Augustus - Version 3.3.3** - http://bioinf.uni-greifswald.de/augustus/  
This tool has three comparisons with the organism models *E. coli* and *S. aureus* and *H. sapiens*.  

**EasyGene - Version 1.2** - http://www.cbs.dtu.dk/services/EasyGene/  
This tool has two comparisons with the organism models *E. coli - K12* and *S. aureus Mu50*.

 **FGENESB - Version '2016'** - http://www.softberry.com/berry.phtml?topic=fgenesb&group=programs&subgroup=gfindb  
This tool has two comparisons with the organism models *E. coli - K12* and *S. aureus MU50*.   

 **GeneMark - Version 2.5** - http://exon.gatech.edu/GeneMark/gm.cgi  
 This tool has two comparisons with the organism models *E. coli - K12 - MG165* and *S. aureus Mu50*.

**GeneMark.hmm  - Version 3.2.5** -  http://exon.gatech.edu/GeneMark/gmhmmp.cgi
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

**MetaGene - Version '2006'** - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1636498/  
Default options were used.

**MetaGeneAnnotator - Version '2008'** - http://metagene.nig.ac.jp/  
Defaults options were used.

**MetaGeneMark - Version '2010'** - http://exon.gatech.edu/meta_gmhmmp.cgi  
GFF was chosen as output type.

**Prodigal - Version 2.62** - https://github.com/hyattpd/Prodigal  
GFF was chosen as output type.

**TransDecoder - Version 5.5.0** - https://github.com/TransDecoder/TransDecoder/wiki  
Defaults options were used.
