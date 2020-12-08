#!/usr/bin/perl

sub revcomp {
  ($seq)=@_;
  $ret=reverse $seq;
  $ret=~tr/ACTGN/TGACN/;
  return $ret;
}

($fa, $maxmm)=@ARGV;
die "maxmm?" unless $maxmm;
open FA, "<$fa" or die "fa?$fa";
while (<FA>) {
  if (/>(\S+)/) { $chr=$1; next; }
  chomp(); die unless $chr;
  $SEQ{$chr}.=uc $_; 
}
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
#chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-A	chr1	100.00	23	0	0	1	23	10608	10630	6.1e-05	45.0
while (<STDIN>) { next if /^#/;
  ($qid, $sid, $pct, $len, $mm, $gap, $qs, $qe, $ss, $se, $eval, $bits)=split();
  die $qid unless $qid=~/(\S+):(\d+)([+-])(\S+)#([ACTG]{20})-(\S)$/;
  ($chr, $pos, $flag, $hex, $sg, $N)=($1,$2,$3,$4,$5,$6);
  die "$_:$qs<$qe?" unless $qs<=$qe;
  $rev=0;
  $rev=1 if ($ss>$se);
  ($ss, $se)=($se, $ss) if $rev; 
  if (!$rev) {
    $ss-=($qs-1); $se+=(23-$qe);
  } else {
    $ss-=(23-$qe); $se+=($qs-1);
  }
  die "$sid not found? $_" unless exists $SEQ{$sid};
  next if $gap;
  # die "22?: $_" unless (abs($se-$ss)==22);
  $sseq=substr($SEQ{$sid}, $ss-1, 23);
  $sseq=revcomp($sseq) if $rev;
 
  $key="$chr:$pos$flag$hex#$sg";
  $pos-=3 if ($flag=~/-/); $pos2=$pos+22;

  $mm=0;
  @q=split(//, $sg); @s=split(//, $sseq);
  for ($i=0; $i<20; $i++) {
    $mm++ if ($q[$i] ne $s[$i]); 
  }
  next unless ($s[21] eq 'G' && $s[22] eq 'G');
  next if $mm>$maxmm;
  # print;
  print "$key\t$chr:$pos-$pos2\t$sg${N}GG\t$mm\t$sid:$ss-$se\t$sseq\n";
}

__END__
# BLASTN 2.9.0+
# Query: chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A
# Database: ../bin/ecoli
# Fields: query id, alignment length, mismatches, evalue, q. start, q. end, subject id, s. start, s. end, subject strand, subject seq
# 25 hits found
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	20	0	3.83e-05	1	20	chr1	5575	5556	minus	TCGCATCAGGCGCTGAATGC
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	19	0	1.51e-04	1	19	chr1	2894848	2894830	minus	TCGCATCAGGCGCTGAATG
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	20	1	0.009	1	20	chr1	66632	66651	plus	TCGCATCAGGCGTTGAATGC
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	20	1	0.009	1	20	chr1	66717	66736	plus	TCGCATCAGGCGTTGAATGC
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	20	1	0.009	1	20	chr1	2118572	2118553	minus	TCGCATCAGGCACTGAATGC

-bash-4.1$ zcat ecoli.m1.blast.gz | awk '$2==23&&$3==0' | head
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-C 23  0   6.21e-07    1   23  chr1    5575    5553    minus   TCGCATCAGGCGCTGAATGCCGG
chr1:5559+74de39771a#TTCAGCGCCTGATGCGACGC-T 23  0   6.21e-07    1   23  chr1    5559    5581    plus    TTCAGCGCCTGATGCGACGCTGG

chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-A	chr1	100.00	23	0	0	1	23	10608	10630	6.1e-05	45.0
chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-C	chr1	95.65	23	1	0	1	23	10608	10630	1.1e-03	41.0
chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-T	chr1	95.65	23	1	0	1	23	10608	10630	1.3e-03	41.0
chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-G	chr1	95.65	23	1	0	1	23	10608	10630	3.4e-04	43.0
chr1:10619-f41ea5ddf7#GCGGCGCGCCTTTGCAACGG-C	chr1	100.00	23	0	0	1	23	10638	10616	3.6e-05	46.0
chr1:10621+775ddf01eb#GTTGCAAAGGCGCGCCGCGC-A	chr1	95.65	23	1	0	1	23	10621	10643	7.6e-04	42.0
chr1:10621+775ddf01eb#GTTGCAAAGGCGCGCCGCGC-C	chr1	100.00	23	0	0	1	23	10621	10643	3.6e-05	46.0
chr1:10621+775ddf01eb#GTTGCAAAGGCGCGCCGCGC-C	chr1	100.00	16	0	0	8	23	181166	181181	4.6e-01	33.0
chr1:10621+775ddf01eb#GTTGCAAAGGCGCGCCGCGC-C	chrY	95.65	23	1	0	1	23	57215899	57215877	2.1e-04	44.0
chr1:10621+775ddf01eb#GTTGCAAAGGCGCGCCGCGC-C	chrX	95.65	23	1	0	1	23	156029379	156029357	2.1e-04	44.0
# BLAT 36x2 [2009/02/26]
# Query: chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-A
# Database: ../bin/hg38.fa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-A	chr1	100.00	23	0	0	1	23	10608	10630	6.1e-05	45.0
# BLAT 36x2 [2009/02/26]
# Query: chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-C
# Database: ../bin/hg38.fa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
chr1:10608+7ad759074#ACGCAACTCCGCCGTTGCAA-C	chr1	95.65	23	1	0	1	23	10608	10630	1.1e-03	41.0
