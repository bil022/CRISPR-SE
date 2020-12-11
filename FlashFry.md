### Scoring function:
#### CRISPR-SE uses much less time in build index and query. CRISPR-SE provides a script to convert the outputs into the format that can be scored using FlashFry.
#### The scripts provided by FlashFry for index, query and scoring:
```
# Inputs:
#   EcoliE84.fa: reference genome
#.  ecoli_input.fa: query
input=ecoli_input

# Build index for reference genome
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 index \
 --tmpLocation ./tmp \
 --database ecoli_cas9ngg_database \
 --reference EcoliE84.fa \
 --enzyme spcas9ngg

# Search off-targets
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 discover \
 --database ecoli_cas9ngg_database \
 --fasta $input.fa \
 --output $input.output

# Scoring
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 score \
 --input $input.output \
 --output $input.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_cas9ngg_database
```

#### How to run FlashFry scoring function with CRISPR-SE
```
# Build index for reference genome, one time only
$BIN/se --index -r EcoliE84
# Build index for user input
$BIN/se --index -r $input
# Create sequence context for input
# Example: 31f78d:CGATGCGGCAGATTGCGCCA	0	ecoli_input	11	30	20M	*	0	0	CGATGCGGCAGATTGCGCCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.794850
./flashfry.pl --context --query $input > $input.ctx

# Convert gRNA into fa format
sed 's/:/ /' $input.ref | awk '{print ">"$1"\n"$2}' > $input.se.fa

# Convert output into FlashFry format
$BIN/se --build -m 5 -v -r $REF -q $input | ./flashfry.pl --format | sort | ./flashfry.pl --parse --ref $REF --query $input.se | ./flashfry.pl --discover --query $input | sort -k1,1 -k2,2n > $input.se.output

# Run FlashFry for scoring function
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 score \
 --input $input.se.output \
 --output $input.se.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_cas9ngg_database
```
