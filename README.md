# CRISPR-SE
## CRISPR-SE is an accurate and fast guide RNA off-target sites search engine for CRISPR design. *The scripts and full datasets are included in [demo](demo/) and [benchmark](benchmark/) folder in [online](http://renlab.sdsc.edu/CRISPR-SE/) and Zenodo(TBD) version.*

## Source code:
<!--
Source: /projects/ps-renlab/bil022/public_html/CREST-web/rev/ref/se_rev/
-->
```
git clone https://github.com/bil022/CRISPR-SE
cd CRISPR-SE
make
```

## Inputs:
```
# demo/ecoli.fa: provided reference genome, can be replaced with mm10, hg38, etc
# demo/crispor_input.fa: provided test input for CRISPOR
# demo/flashfry_input.fa: provided test input for FlashFry
```

## Genome-wide guide RNA design:
1. Build index for reference genome:
```
$ ./se --index -r ecoli
```
>ecoli.h: header files in SAM format<br/>
>ecoli.rep: guide RNA repeats<br/>
>ecoli.ref: unique reference guide RNA<br/>
>ecoli.idx: guide RNA index

2. Genome wide gRNA design:
```
$ ./se --build -r ecoli -p 48
```
>Load 529028 gRNAs in ecoli<br/>
>ecoli: 0.174117 seconds<br/>
>p9: 122.671851 seconds<br/>
>...<p/>
>ecoli.mm: guide RNA that pass off-target filtering criteria

3. Convert to bam format:
```
$ cat ecoli.h ecoli.mm | samtools view -Sb - | samtools sort > ecoli.bam
$ samtools index ecoli.bam
```
>ecoli.bam: guide RNA in BAM format

4. View guide RNAs:
```
$ samtools view ecoli.bam ecoli:1000-1005
```

>#id:gRNA strand chr pos qual 20M * 0 0 seq efficiency_score<br/>
>40b77d8107:GCAACAATCGGCGCGTAAAC	16	ecoli	986	30	20M	*	0	0	GTTTACGCGCCGATTGTTGC	IIIIIIIIIIIIIIIIIIII	EF:f:0.532925<br/>
>a337ae8d77:GCGCCGATTGTTGCGAGATT	0	ecoli	992	30	20M	*	0	0	GCGCCGATTGTTGCGAGATT	IIIIIIIIIIIIIIIIIIII	EF:f:0.043497<br/>
>...

## How to replace CRISPOR search engine:
1. Inteteral processing in CRISPOR:
```
$ bwa aln -o 0 -m 1980000 -n 4 -k 4 -N -l 20 ecoli.fa crispor_input.fa > crispor_input.sa
$ bwa samse -n 60000 ecoli.fa crispor_input.sa crispor_input.fa | ./xa2multi.pl
```

>ecoli:254221	16	ecoli	254221	12	20M	*	0	0	GAAGGCGATTAAACGCCATC	*	XT:A:U	NM:i:0	X0:i:1	X1:i:13	XM:i:0	XO:i:0	XG:i:0	MD:Z:20	XA:Z:ecoli,+4564806,20M,4;ecoli,+3811590,20M,4;ecoli,+2143427,20M,4;ecoli,+2977513,20M,4;ecoli,-1675737,20M,4;ecoli,+3741225,20M,4;ecoli,-805488,20M,4;ecoli,-92861,20M,4;ecoli,-969859,20M,4;ecoli,+782981,20M,4;ecoli,-2349904,20M,4;ecoli,-1156520,20M,4;<br/>
>...<br/>

2. Replace with CRISPR-SE
```
$ ./se --index -sr crispor_input
$ ./se --build -m 5 -v -r ecoli -q crispor_input | ./crispor-se.pl --format | sort | ./crispor-se.pl --parse --ref ecoli --query crispor_input
```

>ecoli:254221	16	ecoli	254221	30	20M	*	0	0	GAAGGCGATTAAACGCCATC	IIIIIIIIIIIIIIIIIIII	EF:f:0.476010	NM:i:0<br/>
>ecoli:254221	0	ecoli	254221	30	20M	*	0	0	GAAGGCGATTAAACGCCATC	IIIIIIIIIIIIIIIIIIII	EF:f:0.454607	NM:i:4<br/>
>ecoli:254221	0	ecoli	782981	30	20M	*	0	0	GATGGCCATGAATGGCCTTC	IIIIIIIIIIIIIIIIIIII	EF:f:0.054218	NM:i:4<br/>

3. The results will further be processed by CRISPOR web tools.

>For example, query GATGGCGTTTAATCGCCTTCCGG at Chromosome:254218-254240 (ecoli) reverse strand. Updated [CRISPOR with CRISPR-SE](http://renlab.sdsc.edu/CRISPR-SE/crispor/crispor.py) reported two off-targets, where origin [CRISPOR](http://crispor.tefor.net/crispor.py) only report the first one due to limits in heuristic algorithm:
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

## Scoring with FlashFry:
<!-- /projects/ps-renlab/bil022/public_html/CRISPR-SE/demo -->
### A. Use FlashFry
1. Build index for reference genome
```
$ java -Xmx4g -jar FlashFry.jar \
 index \
 --tmpLocation ./tmp \
 --database ecoli_db \
 --reference ecoli.fa \
 --enzyme spcas9ngg
```
2. Search off-targets
```
$ java -Xmx4g -jar FlashFry.jar \
 discover \
 --database ecoli_db \
 --fasta flashfry_input.fa \
 --output flashfry_input.output
```

>#contig  start   stop    target  context overflow        orientation     otCount offTargets<br/>
>ecoli   7       30      TGGCGCAATCTGCCGCATCGTGG GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG     OK      RVS     1       TGGCGCAATCTGCCGCATCGTGG_1_0<br/>
>...

3. Scoring
```
$ java -Xmx4g -jar FlashFry.jar \
 score \
 --input flashfry_input.output \
 --output flashfry_input.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_db
```

>contig  start   stop    target  context overflow        orientation     Doench2014OnTarget      Doench2016CFDScore      dangerous_GC    dangerous_polyT dangerous_in_genome     Hsu2013 basesDiffToClosestHit   closestHitCount otCount<br/>
>ecoli   7       30      TGGCGCAATCTGCCGCATCGTGG GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG     OK      RVS     0.1414161053171064      0.0     NONE    NONE    IN_GENOME=1     100.0   UNK     0       1<br/>


### B. Use CRISPR-SE for faster off-targets sites search

1. Build index for user input
```
$ ./se --index -r flashfry_input
```
2. Create sequence context for input
```
$ ./flashfry.pl --context --query flashfry_input > flashfry_input.ctx
```
>177a31f78d:CGATGCGGCAGATTGCGCCA	0	ecoli	11	30	20M	*	0	0	CGATGCGGCAGATTGCGCCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.794850	PM:Z:CGG	TX:Z:AATCCACGATGCGGCAGATTGCGCCACGGCGCTGT<br/>
>...
3. Convert gRNA into fasta format
```
$ sed 's/:/ /' flashfry_input.ref | awk '{print ">"$1"\n"$2}' > flashfry_input.se.fa
```
>\>177a31f78d<br/>
>CGATGCGGCAGATTGCGCCA<br/>
>...

4. Convert output into FlashFry format
```
$ ./se --build -m 5 -v -r ecoli -q flashfry_input | ./flashfry.pl --format | sort | ./flashfry.pl --parse --ref ecoli --query flashfry_input.se | ./flashfry.pl --discover --query flashfry_input | sort -k1,1 -k2,2n > flashfry_input.se.output
```
>contig	start	stop	target	context	overflow	orientation	otCount	offTargets<br/>
>ecoli	7	30	TGGCGCAATCTGCCGCATCGTGG	GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG	OK	RVS	1	TGGCGCAATCTGCCGCATCGTGG_1_0<br/>
>...<br/>

5. Use FlashFry for scoring
```
$ java -Xmx4g -jar FlashFry.jar \
 score \
 --input flashfry_input.se.output \
 --output flashfry_input.se.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_db
```
>contig	start	stop	target	context	overflow	orientation	Doench2014OnTarget	Doench2016CFDScore	dangerous_GC	dangerous_polyT	dangerous_in_genome	Hsu2013	basesDiffToClosestHit	closestHitCount	otCount<br/>
>ecoli	7	30	TGGCGCAATCTGCCGCATCGTGG	GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG	OK	RVS	0.1414161053171064	0.0	NONE	NONE	IN_GENOME=1	100.0	UNK	0	1<br/>
> ...<br/>

### C. Compare the different
```
$ diff flashfry_input.output.scored flashfry_input.se.output.scored
```

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
run `make` in benchmark folder to re-run all test cases. In detail:
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

