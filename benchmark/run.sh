#!/bin/bash

input=input
output=output
blastn=blastn

time=/usr/bin/time
limit=10000

function filt {
  awk '/[0-9a-f]$/{s++}/-A$/{s++}{if('$limit'>0&&s>'$limit')exit;print}'
}

for ref in mm10 hg38; do
  for m in 1 2 3 4; do
    echo $ref m$m begin
    if ! [ -e blast/$ref.m$m.blast.tm ]; then
      $time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o blast/$ref.m$m.blast.tm $blastn -task blastn-short -ungapped -num_threads 8 -db ../bin/$ref -query <(zcat fa/$ref.m$m.fasta.gz|filt) -evalue 1 -outfmt "7 delim=	 qseqid length mismatch evalue qstart qend sseqid sstart send sstrand sseq" 2> blast/$ref.m$m.blast.log | gzip > blast/$ref.m$m.blast.gz & sleep 1
    fi

    if ! [ -e bowtie/$ref.m$m.bowtie.tm ]; then
      v=`expr $m + 1`
      if [ $v -gt 3 ]; then v=3; fi
      $time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o bowtie/$ref.m$m.bowtie.tm bowtie -S -a -v $v -f ../bin/$ref <(zcat fa/$ref.m$m.fasta.gz|filt) 2> bowtie/$ref.m$m.bowtie.log | gzip > bowtie/$ref.m$m.bowtie.gz & sleep 1 
    fi
    
    if ! [ -e bowtie2/$ref.m$m.bowtie2.tm ]; then # -a => -k 20
      echo "../bin/bowtie2 -a -x ../bin/$ref -k 100 -N 1 -f <(zcat fa/$ref.m$m.fasta.gz|filt) -L 12 2> bowtie2/$ref.m$m.bowtie2.log | gzip > bowtie2/$ref.m$m.bowtie2.gz"
      $time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o bowtie2/$ref.m$m.bowtie2.tm bowtie2 -a -x ../bin/$ref -k 100 -N 1 -f <(zcat fa/$ref.m$m.fasta.gz|filt) -L 12 2> bowtie2/$ref.m$m.bowtie2.log | gzip > bowtie2/$ref.m$m.bowtie2.gz & sleep 1
    fi

    if ! [ -e bwa/$ref.m$m.bwa.tm ]; then
      ($time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o bwa/$ref.m$m.bwa.tm bwa aln -a -o 0 -n 4 -k 4 -N ../bin/$ref.fa <(zcat fa/$ref.m$m.fasta.gz|filt) 2> bwa/$ref.m$m.bwa.log > bwa/$ref.m$m.bwa.sai
      bwa samse ../bin/$ref.fa bwa/$ref.m$m.bwa.sai <(zcat fa/$ref.m$m.fasta.gz|filt) 2> bwa/$ref.m$m.bwa.samse.log | samtools view -Sb - > bwa/$ref.m$m.bwa.bam) & sleep 1
    fi

    if ! [ -e blat/$ref.m$m.blat.tm ]; then
      zcat fa/$ref.m$m.fasta.gz | filt > blat/$ref.m$m.fa 
      ($time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o blat/$ref.m$m.blat.tm blat -out=blast9 -stepSize=5 -repMatch=2253 -oneOff=1 -maxGap=0 -minScore=0 -minIdentity=4 ../bin/$ref.fa blat/$ref.m$m.fa blat/$ref.m$m.blat >& blat/$ref.m$m.blat.log
      gzip blat/$ref.m$m.blat
      rm blat/$ref.m$m.fa ) & sleep 1
    fi

    wait

    echo $ref m$m done
  done
  exit
done
