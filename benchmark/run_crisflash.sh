#!/bin/bash

time=/usr/bin/time
limit=10000
function filt {
  awk '/[0-9a-f]$/{s++}/-A$/{s++}{if('$limit'>0&&s>'$limit')exit;print}'
}

BASE=/projects/ps-renlab/bil022/CRISPR/CREST-HPC/refD/rev
CRISFLASH="$BASE/crisflash/crisflash"

for ref in mm10 hg38; do
  for m in 1 2 3 4; do
    echo $ref m$m begin
    m1=`expr $m + 1`
    if ! [ -e crisflash/$ref.m$m.crisflash.tm ]; then
      echo "( $time -f \"%e real,%U user,%S sys,%P CPU,%K mem(K): %C\" -o crisflash/$ref.m$m.crisflash.tm $CRISFLASH -g ref/$ref.fa -s ref/$ref.m$m.fa -o crisflash/$ref.m$m.crisflash -m $m1 -t 8 -C > crisflash/$ref.m$m.crisflash.log; gzip crisflash/$ref.m$m.crisflash ) & sleep 1"
      ( $time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o crisflash/$ref.m$m.crisflash.tm $CRISFLASH -g ref/$ref.fa -s ref/$ref.m$m.fa -o crisflash/$ref.m$m.crisflash -m $m1 -t 8 -C > crisflash/$ref.m$m.crisflash.log; gzip crisflash/$ref.m$m.crisflash ) & sleep 1
    fi
    wait
    echo $ref m$m done
  done
done

# crisflash/crisflash -g ref/mm10.fa -s ref/mm10.m4.fa -o crisflash/mm10.m4.crisflash -m 5 -t 8 -A > crisflash/mm10.m4.crisflash.log 
