#!/bin/bash

for ref in hg38 mm10; do
  for m in 1 2 3 4; do
    echo -n "$ref $m "
    samtools view $ref.m$m.bwa.bam | ./bwa.pl ../../bin/$ref.fa $m | awk '{print $1,$4}' | sort -u | awk 'BEGIN{a=z=m=0}{a++;if($NF==0)z++;if($NF=='$m')m++;}END{print a,z,m,a-z-m}'
  done
done
