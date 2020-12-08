for ref in mm10 hg38; do
  for m in 1 2 3 4; do
    echo "ref=$ref m=$m"
    if ! [ -e $ref.m$m.se.log ]; then
      /usr/bin/time -f "%e real,%U user,%S sys,%P CPU,%K mem(K): %C" -o $ref.m$m.se.tm ./se --build -r $ref -q $ref.m$m -m$m -p 8 >& $ref.m$m.se.log
    fi
  done
done
