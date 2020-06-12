# Pro_Gene - Prokaryote Genome Anotation Comparison 
Comparison pipeline for Prokaryote Protein Coding Gene Predictors 

For Help: python3 ORF_Compare.py -h   
Example: python3 ORF_Compare.py -t Prodigal -i Prodigal_E-coli.gff -g E-coli  
```buildoutcfg
usage: ORF_Compare.py [-h] -t TOOL -i INPUT_TO_ANALYSE -g GENOME_TO_COMPARE

optional arguments:
  -h, --help            show this help message and exit
  -t TOOL, --tool TOOL  Which tool to compare?
  -i INPUT_TO_ANALYSE, --input_to_analyse INPUT_TO_ANALYSE
                        Location of tool output to compare.
  -g GENOME_TO_COMPARE, --genome_to_compare GENOME_TO_COMPARE
                        Which genome to analyse? Genome files have same prefix
                        - .fa and .gff appended

```

Output = "'genome_to_compare'.csv"   
Each prediction tool has a directory and data handling script which converts the different tool output into a dictionary the 
comparator understands.  
New genomes and tools can be added to the analysis but they must follow the same directory and naming structure as originals.