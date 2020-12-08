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

while (<STDIN>) { next if /^#/;
  ($qid, $len, $mm, $eval, $qs, $qe, $sid, $ss, $se, $strand, $seq)=split();
  die $qid unless $qid=~/(\S+):(\d+)([+-])(\S+)#([ACTG]{20})-(\S)$/;
  ($chr, $pos, $flag, $hex, $sg, $N)=($1,$2,$3,$4,$5,$6);
  die unless $qs<$qe;
  $rev=0;
  $rev=1 if ($ss>$se);
  ($ss, $se)=($se, $ss) if $rev; 
  if (!$rev) {
    $ss-=($qs-1); $se+=(23-$qe);
  } else {
    $ss-=(23-$qe); $se+=($qs-1);
  }
  die unless exists $SEQ{$sid};
  die unless (abs($se-$ss)==22);
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
