# ORForise - Prokaryote Genome Anotation Comparison 
Comparison pipeline for Prokaryote Protein Coding Gene Predictors

To run, you need:   
* Input Genome FASTA and corresponding GFF file with predicted ORFs for the genome you want to analyse.  
* A prediction output from one of the compatible tools for the same genome (New tools can be added).   
* Python 3 and the numpy library.

For Help: python3 ORF_Compare.py -h   
Example: python3 ORF_Compare.py -t Prodigal -g E-coli \
Example with Model Parameter: python3 ORF_Compare.py -t GeneMark -p Staph -g Staph  


```python
usage: ORF_Compare.py [-h] -t TOOL [-p PARAMETERS] -g GENOME_TO_COMPARE

optional arguments:
  -h, --help            show this help message and exit
  -t TOOL, --tool TOOL  Which tool to compare?
  -p PARAMETERS, --parameters PARAMETERS
                        Optional parameters for prediction tool.
  -g GENOME_TO_COMPARE, --genome_to_compare GENOME_TO_COMPARE
                        Which genome to analyse? Genome files have same prefix
                        - .fa and .gff appended


```

Output = "'genome_to_compare'.csv"   
Each prediction tool has a directory and data handling script which converts the different tool output into a dictionary the
comparator understands.  
New genomes and tools can be added to the analysis but they must follow the same directory and naming structure as originals.

# Genomes Compared:
* Escherichia coli K-12* - Strain ER3413 - Assembly ASM80076v1    
* Staphylococcus aureus* - Strain 502A - Assembly ASM59796v1  
* Bacillus subtilis* - Strain BEST7003 - Assembly ASM52304v1  
* Mycoplasma genitalium* - Strain G37 - Assembly ASM2732v1  
* Caulobacter crescentus* - Strain CB15 - Assembly ASM690v1  
* Pseudomonas fluorescens* - Strain UK4 - Assembly ASM73042v1  



# Tools Compared:   
There are two Groups of tools - Those which do require a pre-built model and those which do not. They are listed with the non-default options used:
## Model Based Tools:

**Augustus** - http://bioinf.uni-greifswald.de/augustus/  
This tool has three comparisons with the organism models *E. coli* and *S. aureus* and *H. sapiens*.  

**EasyGene** - http://www.cbs.dtu.dk/services/EasyGene/  
This tool has two comparisons with the organism models *E. coli - K12* and *S. aureus Mu50*.

 **FGENESB** - http://www.softberry.com/berry.phtml?topic=fgenesb&group=programs&subgroup=gfindb  
This tool has two comparisons with the organism models *E. coli - K12* and *S. aureus MU50*.   

 **GeneMark** - http://exon.gatech.edu/GeneMark/gm.cgi  
 This tool has two comparisons with the organism models *E. coli - K12 - MG165* and *S. aureus Mu50*.

**GeneMark HMM** -  http://exon.gatech.edu/GeneMark/gmhmmp.cgi
 This tool has two comparisons with the organism models *E. coli - K12 - MG165* and *S. aureus Mu50*.    

## Self-Training/Non-Model Based Tools

**FragGeneScan** - https://omics.informatics.indiana.edu/FragGeneScan/    
The 'complete' genome option was selected and GFF was chosen as output type.    

**GeneMark HA** - http://exon.gatech.edu/GeneMark/heuristic_gmhmmp.cgi  
GFF was chosen as output type.

**GeneMarkS** - http://exon.gatech.edu/GeneMark/genemarks.cgi  
GFF was chosen as output type.  

**GeneMarkS-2** - http://exon.gatech.edu/GeneMark/genemarks2.cgi   
GFF3 was chosen as output type.  

**GLIMMER-3** - http://ccb.jhu.edu/software/glimmer/index.shtml  
Default parameters from manual were used.  

**MetaGene** - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1636498/  
Default options were used.

**MetaGeneAnnotator** - http://metagene.nig.ac.jp/  
Defaults options were used.

**MetaGeneMark** - http://exon.gatech.edu/meta_gmhmmp.cgi  
GFF was chosen as output type.

**Prodigal** - https://github.com/hyattpd/Prodigal  
GFF was chosen as output type.

**TransDecoder** - https://github.com/TransDecoder/TransDecoder/wiki  
Defaults options were used.
