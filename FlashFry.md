## Scoring function:

### There are many scoring functions available. CRISPR-SE provides a script to convert the outputs into the format that can be scored by FlashFry.

#### Example of using FlashFry to build index for reference genome, search off-targets and scoring. 
#### Inputs: reference genome *EcoliE84.fa* and query *ecoli_input.fa*
```
input=ecoli_input

java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 index \
 --tmpLocation ./tmp \
 --database ecoli_cas9ngg_database \
 --reference EcoliE84.fa \
 --enzyme spcas9ngg

java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 discover \
 --database ecoli_cas9ngg_database \
 --fasta $input.fa \
 --output $input.output

java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 score \
 --input $input.output \
 --output $input.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_cas9ngg_database
```


```
$BIN/se --index -r $input
./flashfry.pl --context --query $input > $input.ctx

#31f78d:CGATGCGGCAGATTGCGCCA	0	ecoli_input	11	30	20M	*	0	0	CGATGCGGCAGATTGCGCCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.794850
sed 's/:/ /' $input.ref | awk '{print ">"$1"\n"$2}' > $input.se.fa

$BIN/se --build -m 5 -v -r $REF -q $input | ./flashfry.pl --format | sort | ./flashfry.pl --parse --ref $REF --query $input.se | ./flashfry.pl --discover --query $input | sort -k1,1 -k2,2n > $input.se.output

java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 score \
 --input $input.se.output \
 --output $input.se.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_cas9ngg_database
```
