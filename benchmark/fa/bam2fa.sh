for ref in hg38 mm10; do
  for m in 1 2 3 4; do
    echo $ref m$m
    if ! [ -e $ref.m$m.fa.gz ]; then
      samtools view ../mm/$ref.m$m.srt.bam | sed 's/:/ /' | awk '{s="+";if($3)s="-"; print ">"$4":"$5s$1"\n"$2"NGG"}' | gzip > $ref.m$m.fa.gz &
    fi
    if ! [ -e $ref.m$m.fasta.gz ]; then
      samtools view ../mm/$ref.m$m.srt.bam | sed 's/:/ /' | awk '{s="+";if($3)s="-"; 
        print ">"$4":"$5s$1"#"$2"-A\n"$2"AGG"
        print ">"$4":"$5s$1"#"$2"-C\n"$2"CGG"
        print ">"$4":"$5s$1"#"$2"-T\n"$2"TGG"
        print ">"$4":"$5s$1"#"$2"-G\n"$2"GGG"
      }' | gzip > $ref.m$m.fasta.gz &
    fi
  done
done

wait

echo done
# 7839df1876:TCGCATCAGGCGCTGAATGC	16	chr1	5556	30	20M	=	2118553	16	GCATTCAGCGCCTGATGCGA	TCGCATCAGGCACTGAATGC	EF:f:0.089798
