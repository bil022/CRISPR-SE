#!/usr/bin/perl

sub revcomp {
  ($seq)=@_;
  $ret=reverse $seq;
  $ret=~tr/ACTGN/TGACN/;
  return $ret;
}

($fa, $maxmm)=@ARGV;
die "maxmm?" unless $maxmm;

#open FA, "<$fa" or die "fa?$fa";
#while (<FA>) {
#  if (/>(\S+)/) { $chr=$1; next; }
#  chomp(); die unless $chr;
#  $SEQ{$chr}.=uc $_; 
#}

#contig	start	stop	target	context	overflow	orientation	otCount	offTargets
#chr1:3000239-38f9403c5#CCAGGAAAACCTGGATGAAA-A	0	23	CCTTTTCATCCAGGTTTTCCTGG	NONE	OK	RVS	37	CCATTTCATCCAGGTTTTCCAGG_34_1<chr18:5798614^F|chr15:33042120^R|chr7:42362864^R|chr14:18750824^R|chr3:62819754^F|chr13:33800250^R|chr15:33243715^F|chr16:64088805^F|chr1:83604327^R|chr2:17178344^F|chr14:15893784^F|chr17:4063464^R|chr16:36974082^R|chr12:114242807^F|chrX:108190939^F|chr6:141732013^R|chr1:159113285^R|chr14:98227753^R|chr2:149201012^R|chr15:65651729^F|chr9:55825249^R|chr2:86723807^F|chr6:27186515^F|chr14:113820449^R|chrX:67392699^F|chr9:109252418^F|chr10:6920750^R|chr3:11679960^F|chr8:18030769^R|chr4:78301792^R|chr15:69125952^R|chr7:14373216^R|chr9:6098121^F|chrX:24029122^R>,CCATTTCATCCAGGTTTTCCTGG_1_1<chr1:3000235^F>,CCTTTTCATCCAGGTTCTCCAGG_1_1<chr4:143677276^R>,CCTTTTCATTCAGGTTTTCCAGG_1_1<chr3:106698061^F>

while (<STDIN>) { next if /^contig\s/;
  ($qid, $qs, $qe, $seq, $ctx, $ovf, $strand, $cnt, $lst)=split();
  die $qid unless $qid=~/(\S+):(\d+)([+-])(\S+)#([ACTG]{20})-(\S)$/;
  warn "qs:$qs qe:$qe\t$_" unless ($qs==0 && $qe==23);
  ($chr, $pos, $flag, $hex, $sg, $N)=($1,$2,$3,$4,$5,$6);
  $rev=0; $rev=1 if ($strand=~/RVS/);
  # die "$qid not found? $_" unless exists $SEQ{$qid};
 
  $key="$chr:$pos$flag$hex#$sg";
  $pos-=3 if ($flag=~/-/); $pos2=$pos+22;

  $base=0; $base=3 if $rev;
  while ($lst=~/([ACTG]{23})_\d+_(\d+)/g) {
    ($offtgt, $mmx)=($1, $2); $mm=0;
    @q=split(//, $offtgt); @s=split(//, $seq);
    for ($i=0; $i<20; $i++) {
      $mm++ if ($q[$i+$base] ne $s[$i+$base]); 
    }
    # warn "$mm!=$mmx: [$offtgt:$mmx] $_" unless $mm==$mmx;
    next if $mm!=$mmx;
    # if ($rev==0) { warn "GG\$? $_" unless $offtgt=~/GG$/; } else { warn "^CC? $_" unless $offtgt=~/^CC/; }
    if ($rev==0) { next unless $offtgt=~/GG$/; } else { next unless $offtgt=~/^CC/; }
    # print;
    print "$key\t$chr:$pos-$pos2\t$sg${N}GG\t$mm\t$offtgt\n";
  }
}
