
# Genome editing tool upgrades:

## Replacing CRISPOR search engine (BWA) with CRSIPR-SE:

1. CRISPOR use BWA for mismatch search with following parameters:
```
bwa aln -o 0 -m 1980000 -n 4 -k 4 -N -l 20 EcoliE84.fa $input.fa > $input.sa
bwa samse -n 60000 EcoliE84.fa $input.sa $input.fa 
```
2. Use CRISPR-SE (se) to generate same output format as BWA: 
```
se --index -sr $input
se --build -m 5 -v -r EcoliE84 -q $input | crispor-se.pl --format | sort | crispor-se.pl --parse --ref EcoliE84 --query $input 
```
3. Example:
Query GATGGCGTTTAATCGCCTTCCGG at Chromosome:254218-254240 (Ecoli) reverse strand. Updated [CRISPOR with CRISPR-SE](http://renlab.sdsc.edu/CRISPR-SE/crispor/crispor.py) reported two off-targets, where origin [CRISPOR](http://crispor.tefor.net/crispor.py) only report the first one due to limits in heuristic algorithm:
```
guide 1:    GATGGCGTTTAATCGCCTTC CGG
off-target: GATGGCCATGAATGGCCTTC AGG
                  ** *   *      
CFD Off-target score: 0.000000
MIT Off-target score: 0.13
Position: Chromosome:782981-783003:+
Distance from target: 0.529 Mbp
```
```
guide 2:    GATGGCGTTTAATCGCCTTC CGG
off-target: GAAGGCGATTAAACGCCATC CGG
              *    *    *    *  
CFD Off-target score: 0.263736
MIT Off-target score: 0.12
Position: Chromosome:254221-254243:+
Distance from target: 0.000 Mbp
```
