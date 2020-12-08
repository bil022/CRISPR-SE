# CRISPR-SE benchmark:

The CRISPR-SE benchmark evaluate the accuracies and speeds of existing mismatch search engines for CRISPR design. We benchmarked 8 search methods: BLAST, BLAT, Bowtie, Bowtie2, BWA, FlashFry, CrisFlash and CRISPR-SE. The full dataset is available http://renlab.sdsc.edu/CRISPR-SE/benchmark/ and also Zenodo.

## Inputs:
```
fa/$ref.m$m.fasta.gz
```
* The guide RNAs that have at least $m mismatches (1-4) to another gRNA in reference genome ($ref, ex: mm10, hg38)

## Program:
```
# program to run BLAST, BLAT, Bowtie, Bowtie2 and BWA 
run.sh
crispr-se/run.sh
run_crisflash.sh
run_flashfry.sh
$prog/$prog.sh
```
## Output:
```
$prog/$prog.txt
$prog/$ref.m$m.$prog.tm
```
