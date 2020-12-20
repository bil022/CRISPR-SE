#!/bin/bash
# git clone https://github.com/bil022/CRISPR-SE
# cd CRISPR-SE
# make

# clean:
# for x in ecoli flashfry_input crispor_input; do rm -f $x.h $x.rep $x.ref $x.idx; done
# rm -f ecoli.mm ecoli.bam ecoli_db flashfry_input.ctx flashfry_input.se.output

echo "-- CRISPR-SE --"
echo
echo Build ecoli index
if ! [ -e "ecoli.idx" ]; then
  ./se --index -r ecoli
else
  echo Found ecoli.idx
fi
echo

echo Search ecoli genome-wide gRNA
if ! [ -e "ecoli.mm" ]; then
  ./se --build -r ecoli -p 48
fi
  echo Found ecoli.mm
echo

echo Convert to bam format
if ! [ -e "ecoli.bam" ]; then
  cat ecoli.h ecoli.mm | samtools view -Sb - | samtools sort > ecoli.bam
  samtools index ecoli.bam
else
  echo Found ecoli.bam
fi
echo

echo samtools view ecoli.bam ecoli:1000-1005
samtools view ecoli.bam ecoli:1000-1005
echo

echo "-- CRISPOR --"
echo Run with bwa
echo bwa aln
if ! [ -e crispor_input.sa ]; then
  bwa aln -o 0 -m 1980000 -n 4 -k 4 -N -l 20 ecoli.fa crispor_input.fa > crispor_input.sa
else
  echo Found crispor_input.sa
fi
echo

echo bwa samse
bwa samse -n 60000 ecoli.fa crispor_input.sa crispor_input.fa 2> crispor_samse.log | ./xa2multi.pl
echo
echo The bwa outputs will further be processed/filtered by CRISPOR web tools.
echo

echo Run with CRISPR-SE
if ! [ -e crispor_input.idx ]; then
  ./se --index -sr crispor_input
else
  echo Found crispor_input.idx
fi
echo

echo Search with CRISPR-SE
./se --build -m 5 -v -r ecoli -q crispor_input 2> se_build.log | ./crispor-se.pl --format | sort | ./crispor-se.pl --parse --ref ecoli --query crispor_input
echo 
echo The CRISPR-SE outputs will further be processed/filtered by CRISPOR web tools.
echo

echo -- FlashFry --
if ! [ -e ecoli_db ]; then
java -Xmx4g -jar FlashFry.jar \
 index \
 --tmpLocation ./tmp \
 --database ecoli_db \
 --reference ecoli.fa \
 --enzyme spcas9ngg
else
  echo Found ecoli_db
fi
echo

echo FlashFry discover
java -Xmx4g -jar FlashFry.jar \
 discover \
 --database ecoli_db \
 --fasta flashfry_input.fa \
 --output flashfry_input.output > flashfry_discover.log
echo

echo FlashFry score
java -Xmx4g -jar FlashFry.jar \
 score \
 --input flashfry_input.output \
 --output flashfry_input.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_db > flashfry_score.log
echo

echo flashfry_input.output.scored
cat flashfry_input.output.scored
echo

echo -- CRISPR-SE --
if ! [ -e flashfry_input.idx ]; then
  ./se --index -r flashfry_input
else
  echo Found flashfry_input.idx
fi
echo

echo Create flashfry context file
if ! [ -e flashfry_input.ctx ]; then
 ./flashfry.pl --context --query flashfry_input > flashfry_input.ctx
 sed 's/:/ /' flashfry_input.ref | awk '{print ">"$1"\n"$2}' > flashfry_input.se.fa
else
  echo Found flashfry_input.ctx
fi
echo

echo Search off-targets sites
if ! [ -e flashfry_input.se.output ]; then
  ./se --build -m 5 -v -r ecoli -q flashfry_input | ./flashfry.pl --format | sort | ./flashfry.pl --parse --ref ecoli --query flashfry_input.se | ./flashfry.pl --discover --query flashfry_input | sort -k1,1 -k2,2n > flashfry_input.se.output
echo
  echo Found flashfry_input.se.output
fi
echo

echo Score with CRISPR-SE
java -Xmx4g -jar FlashFry.jar \
 score \
 --input flashfry_input.se.output \
 --output flashfry_input.se.output.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
 --database ecoli_db > flashfry_se_score.log
echo

echo flashfry_input.se.output.scored
cat flashfry_input.se.output.scored
echo

echo Compare flashfry_input.output.scored flashfry_input.se.output.scored
diff flashfry_input.output.scored flashfry_input.se.output.scored
if [ $? -ne 0 ]; then
  echo Different
else
  echo Identical 
fi
echo Done
