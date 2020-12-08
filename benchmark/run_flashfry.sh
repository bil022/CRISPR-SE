#!/bin/bash

time=/usr/bin/time
limit=1000000
function filt {
  awk '/[0-9a-f]$/{s++}/-A$/{s++}{if('$limit'>0&&s>'$limit')exit;print}'
}

BASE=/projects/ps-renlab/bil022/CRISPR/CREST-HPC/refD/tools/flashfry
BASE=flashfry
JAR=$BASE/FlashFry-assembly-1.9.3.jar
JAVA=/usr/lib/jvm/java-1.8.0/bin/java
FLASHFRY="$JAVA -Xmx32g -jar $JAR"


for ref in mm10 hg38; do
  DB=$BASE/${ref}_ref_db
  if ! [ -e $DB ]; then
    echo $DB
    if ! [ -e $BASE/$ref.fa.gz ]; then echo $ref.fa.gz not found; exit; fi
    [ -e $BASE/${ref}_tmp ] || mkdir $BASE/${ref}_tmp
    date
    echo "$FLASHFRY index --tmpLocation $BASE/${ref}_tmp  --database $DB  --reference $BASE/$ref.fa.gz  --enzyme spcas9ngg"
    $FLASHFRY index --tmpLocation $BASE/${ref}_tmp  --database $DB  --reference $BASE/$ref.fa.gz  --enzyme spcas9ngg
    date
  fi
done

for ref in mm10 hg38 ecoli; do
  DB=$BASE/${ref}_ref_db
  for m in 1 2 3 4; do
    echo $ref m$m begin

    if ! [ -e flashfry/$ref.m$m.flashfry.tm ]; then
      echo "$FLASHFRY discover --positionOutput --maxMismatch=$m --database flashfry/${ref}_ref_db --fasta <(zcat ../fasta/$ref.m$m.fasta.gz|filt) --output flashfry/$ref.m$m.flashfry"
      ( $time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o flashfry/$ref.m$m.flashfry.tm $FLASHFRY discover --positionOutput --maxMismatch=$m --database flashfry/${ref}_ref_db --fasta <(zcat ../fasta/$ref.m$m.fasta.gz|filt) --output flashfry/$ref.m$m.flashfry >& flashfry/$ref.m$m.flashfry.log; gzip flashfry/$ref.m$m.flashfry ) & sleep 1
    fi

    wait
    echo $ref m$m done
  done
done
