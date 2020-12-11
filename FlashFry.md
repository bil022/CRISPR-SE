# Scoring function:

## There are many scoring functions available:

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
