# CRISPR-SE benchmark:

The CRISPR-SE benchmark evaluate the accuracies and speeds of existing mismatch search engines for CRISPR design. We benchmarked 8 search methods: BLAST, BLAT, Bowtie, Bowtie2, BWA, FlashFry, CrisFlash and CRISPR-SE. The full dataset is available at [benchmark](http://renlab.sdsc.edu/CRISPR-SE/benchmark/) and also submitted to Zenodo.

## Prerequisites:
The benchmark was performed in a Linux platform, the programs should be installed and excutable as: blastn, blat, bowtie, bowtie2, bwa, FlashFry-assembly-1.9.3.jar(provided), crisflash and se. Otherwise a path should be provided for the executables.

## Inputs:
Datasets with the guide RNAs that have at least $m(1-4) mismatches to another gRNA in reference genome ($ref: mm10, hg38)
```
fa/$ref.m$m.fasta.gz
```

## Scripts:
run `make` to re-run all test cases. In detail:
```
# run BLAST, BLAT, Bowtie, Bowtie2 and BWA 
./run.sh
# run CRISPR-SE
./crispr-se/run.sh
# run CrisFlash
./run_crisflash.sh
# run FlashFry
./run_flashfry.sh

# In each script, we use the following template to run the program:
$time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o $prog/$ref.m$m.$prog.tm $prog <paramemters> <$ref> <$input>

# To calculate the accuracies, #prog includes BLAST, BLAT, Bowtie, Bowtie2, BWA, CRISPR-SE, CrisFlash and FlashFry
$prog/$prog.sh
```
## Outputs:
```
# The total number of used CPU times
$prog/$ref.m$m.$prog.tm
# The accuracies are evaluated by each $prog/$prog.sh 
$prog/$prog.txt
```
## Example:
```
# To test reference genome mm10 with up to 3 mismatch (m=3) using blast, the input is fa/mm10.m3.fasta.gz
>chr1:3002096+ab885bae99#CTCTTGTTGTCCATATGTTT-A
CTCTTGTTGTCCATATGTTTAGG
...

# run.sh: perform blast search
$time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o blast/$ref.m$m.blast.tm $blastn -task blastn-short -ungapped -num_threads 8 -db ../bin/$ref -query <(zcat fa/$ref.m$m.fasta.gz|filt) -evalue 1 -outfmt "7 delim=qseqid length mismatch evalue qstart qend sseqid sstart send sstrand sseq" 2> blast/$ref.m$m.blast.log | gzip > blast/$ref.m$m.blast.gz

# The outputs blast/mm10.m3.blast.tm indicating total time is 15083 seconds
15083.39 real ...

# blast.sh read the blast/mm10.m1.blast.gz and output the accuracies to blast.txt
blast/blast.sh

# blast.txt: mm10 #mismatches #total_outputs #zero_mismatch #3_mismatches #unexpected
mm10 3 13055 10000 3055 0

# blast/blast.txt shows the accuracies of 31%
```
