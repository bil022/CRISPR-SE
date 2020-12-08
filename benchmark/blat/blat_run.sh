#!/bin/bash

for ref in hg38 mm10; do
  for m in 1 2 3 4; do
    zcat ../fa/$ref.m$m.fasta.gz | head -n 80000 > $ref.m$m.fa 
    ./blat -out=blast9 -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 ../../bin/$ref.fa $ref.m$m.fa $ref.m$m.psl
    gzip $ref.m$m.psl
    rm $ref.m$m.fa
  done
done 
