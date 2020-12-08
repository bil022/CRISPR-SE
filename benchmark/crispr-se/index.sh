# zcat fa/mm10.m1.fa.gz | grep NGG | sed 's/NGG//' | head -n 1000000 > mm10.m1.1M.fa
# ../se --index -sr mm10.m1.1M

time=/usr/bin/time
limit=10000

function filt {
  awk '/[0-9a-f]$/{s++}/-A$/{s++}{if('$limit'>0&&s>'$limit')exit;print}'
}

for ref in hg38 mm10; do
  for m in 1 2 3 4; do
    echo $ref m$m begin
    if ! [ -e $ref.m$m.idx.log ]; then
      zcat ../fa/$ref.m$m.fa.gz|filt|grep NGG|sed 's/NGG//' > $ref.m$m.fa
      $time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o $ref.m$m.idx.tm ./se --index -sr $ref.m$m >& $ref.m$m.idx.log
    fi
  done
done

