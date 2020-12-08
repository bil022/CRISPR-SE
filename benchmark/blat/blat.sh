#!/bin/bash

for ref in mm10 hg38; do
  for m in 1 2 3 4; do
    echo -n "$ref $m "
    zcat $ref.m$m.blat.gz | ./blat.pl ../../bin/$ref.fa $m | awk '{print $1,$4}' | sort -u | awk '{a++;if($NF==0)z++;if($NF=='$m')m++;}END{print a,z,m,a-z-m}'
  done
done
