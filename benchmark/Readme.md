# CRISPR-SE benchmark:

The CRISPR-SE benchmark evaluate the accuracies and speeds of existing mismatch search engines for CRISPR design. We benchmarked 8 search methods: BLAST, BLAT, Bowtie, Bowtie2, BWA, FlashFry, CrisFlash and CRISPR-SE. The full dataset is available http://renlab.sdsc.edu/CRISPR-SE/benchmark/ and also Zenodo.

## Inputs:
```
# The guide RNAs that have at least $m mismatches (1-4) to another gRNA in reference genome ($ref, ex: mm10, hg38)
fa/$ref.m$m.fasta.gz
```

## Programs:
```
# Scripts to run BLAST, BLAT, Bowtie, Bowtie2 and BWA 
run.sh
# Scripts to run CRISPR-SE
crispr-se/run.sh
# Scripts to run CrisFlash
run_crisflash.sh
# Scripts to run FlashFry
run_flashfry.sh

# Scripts to calculate the accuracies, #prog includes BLAST, BLAT, Bowtie, Bowtie2, BWA, CRISPR-SE, CrisFlash and FlashFry
$prog/$prog.sh
```
## Outputs:
```
# The total number of used CPU times
$prog/$ref.m$m.$prog.tm
# The accuracies
$prog/$prog.txt
```
