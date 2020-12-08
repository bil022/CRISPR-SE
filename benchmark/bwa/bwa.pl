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

#chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	16	chr1	5553	255	23M	*	0	0	CCTGCATTCAGCGCCTGATGCGA	IIIIIIIIIIIIIIIIIIIIIII	XA:i:1	MD:Z:2G20	NM:i:1
while (<STDIN>) { next if /^@/;
  #($qid, $len, $mm, $eval, $qs, $qe, $sid, $ss, $se, $strand, $seq)=split();
  ($rid, $flag, $chr, $pos, $n, $cigar, $a, $b, $c, $seq)=split();
  die $rid unless $rid=~/(\S+):(\d+)([+-])(\S+)#([ACTG]{20})-(\S)$/;
  ($chr0, $pos0, $flag0, $hex, $sg, $N)=($1,$2,$3,$4,$5,$6);
  die unless ($flag==0 || $flag==16);
  die unless ($cigar=~/23M/); 
  die unless exists $SEQ{$chr};
  $sseq=substr($SEQ{$chr}, $pos-1, 23);
  $sseq=revcomp($sseq) if $flag;
 
  $key="$chr0:$pos0$flag0$hex#$sg";
  $pos0-=3 if ($flag0=~/-/); $pos2=$pos0+22;

  $mm=0;
  @q=split(//, $sg); @s=split(//, $sseq);
  for ($i=0; $i<20; $i++) {
    $mm++ if ($q[$i] ne $s[$i]); 
  }
  $pos22=$pos+22;
  #print;
  #print "$key\t$chr:$pos0-$pos2\t$sg${N}GG\t$mm\t$chr:$pos-$pos22\t$sseq\n";
  next unless ($s[21] eq 'G' && $s[22] eq 'G');
  # next if $mm>$maxmm;
  print "$key\t$chr:$pos0-$pos2\t$sg${N}GG\t$mm\t$chr:$pos-$pos22\t$sseq\n";
}
print STDERR $X;

__END__
@HD	VN:1.0	SO:unsorted
@SQ	SN:chr1	LN:4641652
@PG	ID:Bowtie	VN:0.12.7	CL:"bowtie -S -a -v 2 -f ../bin/ecoli /dev/fd/63"
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	16	chr1	5553	255	23M	*	0	0	CCTGCATTCAGCGCCTGATGCGA	IIIIIIIIIIIIIIIIIIIIIII	XA:i:1	MD:Z:2G20	NM:i:1
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	16	chr1	2894826	255	23M	*	0	0	CCTGCATTCAGCGCCTGATGCGA	IIIIIIIIIIIIIIIIIIIIIII	XA:i:2	MD:Z:2G0A19	NM:i:2
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	0	chr1	66717	255	23M	*	0	0	TCGCATCAGGCGCTGAATGCAGG	IIIIIIIIIIIIIIIIIIIIIII	XA:i:2	MD:Z:12T7C2	NM:i:2
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	0	chr1	66632	255	23M	*	0	0	TCGCATCAGGCGCTGAATGCAGG	IIIIIIIIIIIIIIIIIIIIIII	XA:i:2	MD:Z:12T7C2	NM:i:2
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-A	16	chr1	2118550	255	23M	*	0	0	CCTGCATTCAGCGCCTGATGCGA	IIIIIIIIIIIIIIIIIIIIIII	XA:i:2	MD:Z:2G8T11	NM:i:2
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-C	16	chr1	5553	255	23M	*	0	0	CCGGCATTCAGCGCCTGATGCGA	IIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:23	NM:i:0
chr1:5556-7839df1876#TCGCATCAGGCGCTGAATGC-C	16	chr1	2894826	255	23M	*	0	0	CCGGCATTCAGCGCCTGATGCGA	IIIIIIIIIIIIIIIIIIIIIII	XA:i:1	MD:Z:3A19	NM:i:1
