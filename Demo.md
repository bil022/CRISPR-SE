## Compile CRISPR-SE:
```
$ git clone https://github.com/bil022/CRISPR-SE
$ cd CRISPR-SE
$ make
```
<!--
Source: /projects/ps-renlab/bil022/public_html/CREST-web/rev/ref/se_rev/
-->

## Genome-wide gRNA design:
1. Build index for reference genome:
```
$ ./se --index -r ecoli
```
ecoli.h: header files in SAM format<br/>
ecoli.rep: guide RNA repeats<br/>
ecoli.ref: unique reference guide RNA<br/>
ecoli.idx: guide RNA index<br/>

2. Genome wide gRNA design:
```
$ ./se --build -r ecoli -p 48
```
```
Load 529028 gRNAs in ecoli
ecoli: 0.174117 seconds
p9: 122.671851 seconds
...
```
ecoli.mm: guide RNA that pass off-target filtering criteria

3. Convert to bam format:
```
$ cat ecoli.h ecoli.mm | samtools view -Sb - | samtools sort > ecoli.bam
$ samtools index ecoli.bam
```
ecoli.bam: guide RNA in BAM format

4. View guide RNAs:
```
$ samtools view ecoli.bam ecoli:1000-1005
```
```
#id:gRNA strand chr pos qual 20M * 0 0 seq efficiency_score
40b77d8107:GCAACAATCGGCGCGTAAAC	16	ecoli	986	30	20M	*	0	0	GTTTACGCGCCGATTGTTGC	IIIIIIIIIIIIIIIIIIII	EF:f:0.532925
a337ae8d77:GCGCCGATTGTTGCGAGATT	0	ecoli	992	30	20M	*	0	0	GCGCCGATTGTTGCGAGATT	IIIIIIIIIIIIIIIIIIII	EF:f:0.043497
```

## How to replacing CRISPOR search engine:
```
# 1. Internal processing in CRISPOR:
$ bwa aln -o 0 -m 1980000 -n 4 -k 4 -N -l 20 ecoli.fa crispor_input.fa > crispor_input.sa
$ bwa samse -n 60000 ecoli.fa crispor_input.sa crispor_input.fa | ./xa2multi.pl
```
```
ecoli:254221	16	ecoli	254221	12	20M	*	0	0	GAAGGCGATTAAACGCCATC	*	XT:A:U	NM:i:0	X0:i:1	X1:i:13	XM:i:0	XO:i:0	XG:i:0	MD:Z:20	XA:Z:ecoli,+4564806,20M,4;ecoli,+3811590,20M,4;ecoli,+2143427,20M,4;ecoli,+2977513,20M,4;ecoli,-1675737,20M,4;ecoli,+3741225,20M,4;ecoli,-805488,20M,4;ecoli,-92861,20M,4;ecoli,-969859,20M,4;ecoli,+782981,20M,4;ecoli,-2349904,20M,4;ecoli,-1156520,20M,4;<br/>
```

```
# 2. Replace with CRISPR-SE
$ ./se --index -sr crispor_input
$ ./se --build -m 5 -v -r ecoli -q crispor_input | ./crispor-se.pl --format | sort | ./crispor-se.pl --parse --ref ecoli --query crispor_input
```
```
ecoli:254221	16	ecoli	254221	30	20M	*	0	0	GAAGGCGATTAAACGCCATC	IIIIIIIIIIIIIIIIIIII	EF:f:0.476010	NM:i:0
ecoli:254221	0	ecoli	254221	30	20M	*	0	0	GAAGGCGATTAAACGCCATC	IIIIIIIIIIIIIIIIIIII	EF:f:0.454607	NM:i:4
ecoli:254221	0	ecoli	782981	30	20M	*	0	0	GATGGCCATGAATGGCCTTC	IIIIIIIIIIIIIIIIIIII	EF:f:0.054218	NM:i:4
```
```
# 3. The results will further be processed by CRISPOR web tools.
```
## Scoring with FlashFry:
<!-- /projects/ps-renlab/bil022/public_html/CRISPR-SE/demo -->
### A. Use FlashFry
```
# 1. Build index for reference genome
$ java -Xmx4g -jar FlashFry.jar \
 index \
 --tmpLocation ./tmp \
 --database ecoli_db \
 --reference ecoli.fa \
 --enzyme spcas9ngg
```
```
# 2. Search off-targets
$ java -Xmx4g -jar FlashFry.jar \
 discover \
 --database ecoli_db \
 --fasta flashfry_input.fa \
 --output flashfry_input.output
```
```
contig  start   stop    target  context overflow        orientation     otCount offTargets
ecoli   7       30      TGGCGCAATCTGCCGCATCGTGG GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG     OK      RVS     1       TGGCGCAATCTGCCGCATCGTGG_1_0
...
```
```
# 3. Scoring
$ java -Xmx4g -jar FlashFry.jar \
 score \
 --input flashfry_input.output \
 --output flashfry_input.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_db
```
```
contig  start   stop    target  context overflow        orientation     Doench2014OnTarget      Doench2016CFDScore      dangerous_GC    dangerous_polyT dangerous_in_genome     Hsu2013 basesDiffToClosestHit   closestHitCount otCount
ecoli   7       30      TGGCGCAATCTGCCGCATCGTGG GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG     OK      RVS     0.1414161053171064      0.0     NONE    NONE    IN_GENOME=1     100.0   UNK     0       1
...
```
### B. Use CRISPR-SE
```
# 1. Build index for user input
$ ./se --index -r flashfry_input
```
```
# 2. Create sequence context for input
$ ./flashfry.pl --context --query flashfry_input > flashfry_input.ctx
```
```
177a31f78d:CGATGCGGCAGATTGCGCCA	0	ecoli	11	30	20M	*	0	0	CGATGCGGCAGATTGCGCCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.794850	PM:Z:CGG	TX:Z:AATCCACGATGCGGCAGATTGCGCCACGGCGCTGT
...
```

```
# 3. Convert gRNA into fasta format
$ sed 's/:/ /' flashfry_input.ref | awk '{print ">"$1"\n"$2}' > flashfry_input.se.fa
```
```
>177a31f78d
CGATGCGGCAGATTGCGCCA
...
```

```
# 4. Convert output into FlashFry format
$ ./se --build -m 5 -v -r ecoli -q flashfry_input | ./flashfry.pl --format | sort | ./flashfry.pl --parse --ref ecoli --query flashfry_input.se | ./flashfry.pl --discover --query flashfry_input | sort -k1,1 -k2,2n > flashfry_input.se.output
```
```
contig	start	stop	target	context	overflow	orientation	otCount	offTargets
ecoli	7	30	TGGCGCAATCTGCCGCATCGTGG	GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG	OK	RVS	1	TGGCGCAATCTGCCGCATCGTGG_1_0
...
```
```
# 5. Run FlashFry for scoring
$ java -Xmx4g -jar FlashFry.jar \
 score \
 --input flashfry_input.se.output \
 --output flashfry_input.se.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_db
```
```
 contig	start	stop	target	context	overflow	orientation	Doench2014OnTarget	Doench2016CFDScore	dangerous_GC	dangerous_polyT	dangerous_in_genome	Hsu2013	basesDiffToClosestHit	closestHitCount	otCount
 ecoli	7	30	TGGCGCAATCTGCCGCATCGTGG	GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG	OK	RVS	0.1414161053171064	0.0	NONE	NONE	IN_GENOME=1	100.0	UNK	0	1
 ...
```
### C. Compare the different
```
$ diff flashfry_input.output.scored flashfry_input.se.output.scored
```
```
# Same
```
<!--
Apps:

* se
1. se --index -r mm10
2. se --build -r mm10

* se benchmark
* CRISPOR-SE
* CRISPOR-SE-web
* FlashFry

Inputs/Outputs
-->
