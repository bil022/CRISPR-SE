#!/usr/bin/perl
use Getopt::Long;

$format=0;
$parse=0;
$ref="sacCer3";
$query="hhPawHaI5huGeIlR4Acg";
$context=0;
$discover=0;

GetOptions (
  'format'=>\$format,
  'parse'=>\$parse,
  'ref=s'=>\$ref,
  'query=s'=>\$query,
  'discover'=>\$discover,
  'context'=>\$context
);

if ($format) {
  formatIt();
} elsif ($parse) {
  parse();
} elsif ($context) {
  context();
} elsif ($discover) {
  discover();
} else {
  help();
}

sub help() {
 #print "$0 --format --query $query | sort | $0 --parse --ref $ref\n";
 print "$0 --context --query $query\n";
 print "$0 --format | sort | $0 --parse --ref $ref --query $query\n";
}

sub discover {
=input
query.ctx all gRNAs, may repeat
73a3ed875e:TGCCGCATCGTGGATTGAGC	16	ecoli_input	71	30	20M	*	0	0	GCTCAATCCACGATGCGGCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.035443	PM:Z:TGG	TX:Z:GCAATCTGCCGCATCGTGGATTGAGCTGGGCAAAA
9ce8fb61d7:GCCGCATCGTGGATTGAGCT	16	ecoli_input	70	30	20M	*	0	0	AGCTCAATCCACGATGCGGC	IIIIIIIIIIIIIIIIIIII	EF:f:0.300179	PM:Z:GGG	TX:Z:CAATCTGCCGCATCGTGGATTGAGCTGGGCAAAAA

stdin: all unique gRNAs, include off-targets with NM:i:n>0
7846cdf96e	0	Chromosome	69547	30	20M	*	0	0	TGTCCTGGCGAGTCACATGC	IIIIIIIIIIIIIIIIIIII	EF:f:0.338701	NM:i:0	RP:i:1
7846cdf96e	0	Chromosome	4503200	30	20M	*	0	0	TGTCCAGGCGATTCACAATC	IIIIIIIIIIIIIIIIIIII	EF:f:0.284846	NM:i:4	RP:i:1
7846cdf96e	0	Chromosome	3238012	30	20M	*	0	0	TGTCCTGGCGAGCGAGATTC	IIIIIIIIIIIIIIIIIIII	EF:f:0.059586	NM:i:4	RP:i:1

=output: all gRNAs, append off-targets, 
c		c	c	c			c					-	c	otCount	offTargets
ecoli_input	27	50	GACTCGCCAGGACAGCGCCGTGG	GCATGTGACTCGCCAGGACAGCGCCGTGGCGCAAT	OK	RVS	1	GACTCGCCAGGACAGCGCCGTGG_1_0
ecoli_input	36	59	TGTCCTGGCGAGTCACATGCAGG	CGGCGCTGTCCTGGCGAGTCACATGCAGGATTTTT	OK	FWD	3	TGTCCAGGCGATTCACAATCTGG_1_4,TGTCCTGGCGAGCGAGATTCAGG_1_4,TGTCCTGGCGAGTCACATGCAGG_1_0
=cut
  print "contig\tstart\tstop\ttarget\tcontext\toverflow\torientation\totCount\toffTargets\n";
  while (<>) { if (/^@/) { next; }
    chomp();
    ($id, $strand, $chr, $pos, $qual, $cigar, $a, $b, $c, $seq, $Is, $ef, $nm, $rp)=split(); die $_ unless $rp=~/RP:i:(\d+)/; $rep=$1;
    $seq=revcmp($seq) if $strand&16;
    if (/NM:i:0/) {
      $LN{$id}="${seq}_${rep}_0";
      # warn "$LN{$id}\n";
    }
    die $_ unless $nm=~/NM:i:(\d+)/; $mm=$1;
    if (/NM:i:[1-9]/) {
      die $id unless exists $LN{$id};
      $LN{$id}="${seq}_${rep}_${mm},$LN{$id}";
      # warn "here:$LN{$id} $_\n";
    }
#7846cdf96e	0	Chromosome	69547	30	20M	*	0	0	TGTCCTGGCGAGTCACATGC	IIIIIIIIIIIIIIIIIIII	EF:f:0.338701	NM:i:0	RP:i:1
#7846cdf96e	0	Chromosome	4503200	30	20M	*	0	0	TGTCCAGGCGATTCACAATC	IIIIIIIIIIIIIIIIIIII	EF:f:0.284846	NM:i:4	RP:i:1
#7846cdf96e	0	Chromosome	3238012	30	20M	*	0	0	TGTCCTGGCGAGCGAGATTC	IIIIIIIIIIIIIIIIIIII	EF:f:0.059586	NM:i:4	RP:i:1
  }
  open CTX, "<$query.ctx" or die "$query.ctx";
  while (<CTX>) { chomp();
    ($key, $strand, $qchr, $pos, $qual, $cigar, $a, $b, $c, $seq, $Is, $ef, $pam, $ctx)=split();
    die $_ unless $ctx=~/TX:Z:(\S+)/; $ctx=$1;
    die unless $key=~/(\S+):(\S+)/; ($id, $gRNA)=($1,$2);
    die unless $pam=~/PM:Z:(\S+)/; $pam=$1;
    if (exists $LN{$id}) {
      $s=$pos-1; $s=$pos-4 if ($strand&16);
      $e=$s+23;
      $dir="FWD"; $dir="RVS" if $strand&16;
      $ot=$LN{$id}; $otn=0; $ots="";
      #warn "OT:$ot\n";
      while ($ot=~/([ACTG]+)(_\d+_\d+)/g) {
        ($g, $tag)=($1, $2);
        $otn++; $ots.="," if $ots;
        $ots="$ots$g$pam$tag";
      }
      print "$qchr\t$s\t$e\t$gRNA$pam\t$ctx\tOK\t$dir\t$otn\t$ots\n";
    } else {
      #$warn = sprintf("$id not found $_");
      #warn $warn;
    }
  }
}

# format query for sorting
sub formatIt() {
  #die "query? $query" unless ( -e "$query.mmv");
  #open MMV, "<$query.mmv" or die "$query.mmv";
  #while (<MMV>) {
  while (<>) {
    die $_ unless /^\[\d+\]([0-9A-F]{5}):\d?([0-9A-F]{5})/;
    ($seed, $distal)=($1, $2);
    $key=(hex $seed)*0x100000+(hex $distal);
    print $key; #.":s[$seed]:d[$distal]";
    if (/-([0-9A-F]{5}):\d?([0-9A-F]{5})=(\d+):(\d+)/) {
      ($seed2, $distal2, $mms, $mmd)=($1,$2,$3,$4);
      $key2=(hex $seed2)*0x100000+(hex $distal2);
      $mm=$mms+$mmd;
      print "\t$key2\t$mm";
    }
    print "\n";
  }
}

sub parse() {
  open HEAD, "$ref.h" or die "$ref.h";
  while (<HEAD>) { print; }

  while (<>) { @data=split();
    if (@data==1) { ($key)=@data;
      $KEYS{$key}++; $ANY{$key}++; push(@LST, $key);
    } elsif (@data==3) { ($key, $key2, $mm)=@data;
      if (exists $KEYS{$key}) {
        push(@{$OFF{$key}}, "$key2:$mm");
        # push(@{$MM{$key}}, $mm);
        $ANY{$key2}++;
      } else {
        #$warn = sprintf("%x not found", $key);
        #warn $warn;
      }
    } else {
      die $_;
    }
  }

  # open $ref.ref
  open REF, "$ref.ref" or die "$ref.ref";
  while (<REF>) { chomp();
    die $_ unless /^([0-9a-f]+):/;
    $key=hex($1);
    next unless exists $ANY{$key};
    # save line info if exists $ANY{$key}
    $LN{$key}=$_;
    ($id, $strand, $chr, $pos, $qual)=split();
    # require rep if 10
    $REP{$key}++ if $qual==10;
  }
  foreach $key (keys %ANY) { die $key unless exists $LN{$key}; }
  
  # if require rep
  if (scalar %REP) {
    open RP, "$ref.rep" or die "$ref.rep";
    while (<RP>) {
      die $_ unless /^([0-9a-f]+):/;
      $key=hex($1);
      next unless exists $REP{$key};
      # save REP
      push(@{$REPS{$key}}, $_);
      $NREP{$key}++;
    }
  }

  # update rid from $query.fa
  open FA, "<$query.fa" or die "$query.fa";
  while (<FA>) {
    $rid=$1 if /^>(\S+)/;
    $RID{$1}=$rid if /^([ACTG]+)/;
  }

  foreach $q (@LST) {
    die $q unless exists $LN{$q};
    # how about $q is repeats?
    $RP="\tRP:i:1"; $RP="\tRP:i:$NREP{$q}" if exists $NREP{$q};
    $ln=$LN{$q};
    die $ln unless $ln=~/:(\S+)/; $seq=$1;
    die "$seq:$ln" unless exists $RID{$seq};
    $rid=$RID{$seq};
    @flds=split(/\s+/, $ln);
    shift @flds;
    print "$rid\t".join("\t", @flds)."\tNM:i:0$RP\n";
    foreach $offmm (@{$OFF{$q}}) {
      die $offmm unless $offmm =~/(\S+):(\d+)/;
      ($off, $mm)=($1,$2);
      die "$offmm:$off:$mm" unless exists $LN{$off};
      $RP="\tRP:i:1"; $RP.="\tRP:i:$NREP{$off}" if exists $NREP{$off};
      
      @offln=($LN{$off});
      @offln=@{$REPS{$off}} if exists $REPS{$off};
      foreach $ln (@offln) {
        @flds=split(/\s+/, $ln);
        shift @flds;
        print "$rid\t".join("\t", @flds)."\tNM:i:$mm$RP\n";
      }
    }
  }
}

sub revcmp {
  ($orig)=@_;
  return $orig if $orig=~/^NA$/;
  my $revcomp = reverse $orig;
  $revcomp =~ tr/ATGCatgc/TACGtacg/;
  return $revcomp;
}

sub context {
  open FA, "$query.fa" or die "$query.fa";
  while (<FA>) { chomp(); 
    if (/>(\S+)/) {
      $rid=$1;
    } else {
      die "No seq ID: $_" unless $rid;
      $SEQ{$rid}.=uc $_;
    }
  }
  procContext("$query.ref", 1);
  procContext("$query.rep", 0);
}

sub procContext {
  ($file, $skipRep)=@_;
  open REF, "$file" or die "$file";
  while (<REF>) { chomp();
    ($rid, $strand, $rid, $cut, $score, $cigar, $a, $b, $c, $seq, $Is, $EF)=split(); die $_ unless $EF=~/EF:f:/; 
    next if $skipRep && $score < 30;
    die $rid unless exists $SEQ{$rid};
    $ref=$SEQ{$rid};
    $s=$cut-1;
    if ($strand&16) { #RVS
      $ctx=$s-9;
      $pam=substr($ref, $ctx+6,3);
      $pam=revcmp($pam);
    } else { #FWD
      $ctx=$s-6;
      $pam=substr($ref, $ctx+26,3);
    }
    #print "ctx:$ctx\n";
    $TX="NA";
    if ($ctx>=0) {
      $TX=substr($ref, $ctx, 35);
      $TX="NA" if length $TX<35;
    }
    $TX=revcmp($TX) if $strand&16;
    print "$_\tPM:Z:$pam\tTX:Z:$TX\n";
  }
}

__END__
177a31f78d	0	Chromosome	69521	30	20M	*	0	0	CGATGCGGCAGATTGCGCCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.794850	NM:i:0
5b9df45de8	0	Chromosome	69532	30	20M	*	0	0	ATTGCGCCACGGCGCTGTCC	IIIIIIIIIIIIIIIIIIII	EF:f:0.010988	NM:i:0
5d93b87960	16	Chromosome	69553	30	20M	*	0	0	GGCGAGTCACATGCAGGATT	IIIIIIIIIIIIIIIIIIII	EF:f:0.149034	NM:i:0
7846cdf96e	0	Chromosome	69547	30	20M	*	0	0	TGTCCTGGCGAGTCACATGC	IIIIIIIIIIIIIIIIIIII	EF:f:0.338701	NM:i:0
7846cdf96e	0	Chromosome	4503200	30	20M	*	0	0	TGTCCAGGCGATTCACAATC	IIIIIIIIIIIIIIIIIIII	EF:f:0.284846	NM:i:4
7846cdf96e	0	Chromosome	3238012	30	20M	*	0	0	TGTCCTGGCGAGCGAGATTC	IIIIIIIIIIIIIIIIIIII	EF:f:0.059586	NM:i:4
d7713c5d93	16	Chromosome	69541	30	20M	*	0	0	CGGCGCTGTCCTGGCGAGTC	IIIIIIIIIIIIIIIIIIII	EF:f:0.669980	NM:i:0
d875e6077e	16	Chromosome	69521	30	20M	*	0	0	CGATGCGGCAGATTGCGCCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.141416	NM:i:0
expected
contig	start	stop	target	context	overflow	orientation	otCount	offTargets
ecoli_input	7	30	TGGCGCAATCTGCCGCATCGTGG	GCGCCGTGGCGCAATCTGCCGCATCGTGGATTGAG	OK	RVS	1	TGGCGCAATCTGCCGCATCGTGG_1_0
ecoli_input	10	33	CGATGCGGCAGATTGCGCCACGG	AATCCACGATGCGGCAGATTGCGCCACGGCGCTGT	OK	FWD	1	CGATGCGGCAGATTGCGCCACGG_1_0
ecoli_input	21	44	ATTGCGCCACGGCGCTGTCCTGG	CGGCAGATTGCGCCACGGCGCTGTCCTGGCGAGTC	OK	FWD	1	ATTGCGCCACGGCGCTGTCCTGG_1_0
ecoli_input	27	50	GACTCGCCAGGACAGCGCCGTGG	GCATGTGACTCGCCAGGACAGCGCCGTGGCGCAAT	OK	RVS	1	GACTCGCCAGGACAGCGCCGTGG_1_0
ecoli_input	36	59	TGTCCTGGCGAGTCACATGCAGG	CGGCGCTGTCCTGGCGAGTCACATGCAGGATTTTT	OK	FWD	3	TGTCCAGGCGATTCACAATCTGG_1_4,TGTCCTGGCGAGCGAGATTCAGG_1_4,TGTCCTGGCGAGTCACATGCAGG_1_0
ecoli_input	39	62	AATCCTGCATGTGACTCGCCAGG	GGCAAAAATCCTGCATGTGACTCGCCAGGACAGCG	OK	RVS	1	AATCCTGCATGTGACTCGCCAGG_1_0

__END__
head hhPawHaI5huGeIlR4Acg.fa 
>s0-
GATATAAAAATATATTATAT
>s35+
TTTATATCATTATTTAGACA
>s61-
AAAGCCCTGAAAGCTCACAA
>s76+
AGTCACCATTGTGAGCTTTC
>s77+
GTCACCATTGTGAGCTTTCA
[bli@renlab se]$ head sacCer3.ref hhPawHaI5huGeIlR4Acg.mm hhPawHaI5huGeIlR4Acg.mmv; tail sacCer3.ref sacCer3.rep 
==> sacCer3.ref <==
0:AAAAAAAAAAAAAAAAAAAA	0	chrI	101285	10	20M	*	0	0	AAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIII	EF:f:0.055472
1:CAAAAAAAAAAAAAAAAAAA	0	chrII	316949	10	20M	*	0	0	CAAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIII	EF:f:0.089703
3:GAAAAAAAAAAAAAAAAAAA	16	chrIII	101877	10	20M	*	0	0	TTTTTTTTTTTTTTTTTTTC	IIIIIIIIIIIIIIIIIIII	EF:f:0.217818
5:CCAAAAAAAAAAAAAAAAAA	16	chrXV	675326	30	20M	*	0	0	TTTTTTTTTTTTTTTTTTGG	IIIIIIIIIIIIIIIIIIII	EF:f:0.149894
8:ATAAAAAAAAAAAAAAAAAA	0	chrX	525781	30	20M	*	0	0	ATAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIII	EF:f:0.139083
9:CTAAAAAAAAAAAAAAAAAA	0	chrVIII	551716	30	20M	*	0	0	CTAAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIII	EF:f:0.110271
c:AGAAAAAAAAAAAAAAAAAA	16	chrIII	101878	30	20M	*	0	0	TTTTTTTTTTTTTTTTTTCT	IIIIIIIIIIIIIIIIIIII	EF:f:0.080567
e:TGAAAAAAAAAAAAAAAAAA	16	chrXV	444499	30	20M	*	0	0	TTTTTTTTTTTTTTTTTTCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.132460
30:AAGAAAAAAAAAAAAAAAAA	0	chrII	698073	30	20M	*	0	0	AAGAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIII	EF:f:0.417520
38:ATGAAAAAAAAAAAAAAAAA	0	chrXI	196209	30	20M	*	0	0	ATGAAAAAAAAAAAAAAAAA	IIIIIIIIIIIIIIIIIIII	EF:f:0.224423

==> hhPawHaI5huGeIlR4Acg.mm <==
31d75c1e13:GACATGCAAGCCGCCGCAGA	0	s205-	21	30	20M	*	0	0	GACATGCAAGCCGCCGCAGA	IIIIIIIIIIIIIIIIIIII	EF:f:-1.000000
46c1a29618:ATCATCCTTATTCAAGTCAC	0	s152+	21	30	20M	*	0	0	ATCATCCTTATTCAAGTCAC	IIIIIIIIIIIIIIIIIIII	EF:f:-1.000000
69b8fa71b8:ATGTCAGCTTGGATGTCTTC	0	s243+	21	30	20M	*	0	0	ATGTCAGCTTGGATGTCTTC	IIIIIIIIIIIIIIIIIIII	EF:f:-1.000000
9a6e3e9c6e:TGTCAGCTTGGATGTCTTCT	0	s244+	21	30	20M	*	0	0	TGTCAGCTTGGATGTCTTCT	IIIIIIIIIIIIIIIIIIII	EF:f:-1.000000
9f8ea7f1a9:CTTTCAGGGCTTTGATGGCT	0	s91+	21	30	20M	*	0	0	CTTTCAGGGCTTTGATGGCT	IIIIIIIIIIIIIIIIIIII	EF:f:-1.000000
f79a5e44eb:GTTGACACTGCCTTCTGCGG	0	s215+	21	30	20M	*	0	0	GTTGACACTGCCTTCTGCGG	IIIIIIIIIIIIIIIIIIII	EF:f:-1.000000

==> hhPawHaI5huGeIlR4Acg.mmv <==
[1]132A2:8622A-12222:8628A=2:2
[1]132A2:28622A-13262:86C6A=1:3
[1]132A2:28622A
[1]132A2:28622A-172A2:84026=1:3
[1]132A2:28622A-322A2:0A22A=2:2
[1]20E93:BD78E
[1]20E93:BD78E-23EA3:3E78E=2:2
[1]31D75:C1E13
[1]46C1A:29618
[1]6A73B:A146C
==> sacCer3.ref <==
fffffc816b:GTTCCAATAGGGGGGGGGGG	16	chrXV	94819	30	20M	*	0	0	CCCCCCCCCCCTATTGGAAC	IIIIIIIIIIIIIIIIIIII	EF:f:0.261222
fffffde3ac:AGTTGATGCGGGGGGGGGGG	16	chrX	639717	30	20M	*	0	0	CCCCCCCCCCCGCATCAACT	IIIIIIIIIIIIIIIIIIII	EF:f:0.038218
fffffe80b4:ACGTAAATTGGGGGGGGGGG	16	chrIX	138433	30	20M	*	0	0	CCCCCCCCCCCAATTTACGT	IIIIIIIIIIIIIIIIIIII	EF:f:0.062542
ffffff003f:GGGAAAAAGGGGGGGGGGGG	0	chrV	545179	30	20M	*	0	0	GGGAAAAAGGGGGGGGGGGG	IIIIIIIIIIIIIIIIIIII	EF:f:0.020411
ffffff78eb:GTTGATGCGGGGGGGGGGGG	16	chrX	639716	30	20M	*	0	0	CCCCCCCCCCCCGCATCAAC	IIIIIIIIIIIIIIIIIIII	EF:f:0.032926
ffffffc00f:GGAAAAAGGGGGGGGGGGGG	0	chrV	545180	30	20M	*	0	0	GGAAAAAGGGGGGGGGGGGG	IIIIIIIIIIIIIIIIIIII	EF:f:0.019801
ffffffde3a:TTGATGCGGGGGGGGGGGGG	16	chrX	639715	30	20M	*	0	0	CCCCCCCCCCCCCGCATCAA	IIIIIIIIIIIIIIIIIIII	EF:f:0.002395
fffffff78e:TGATGCGGGGGGGGGGGGGG	16	chrX	639714	30	20M	*	0	0	CCCCCCCCCCCCCCGCATCA	IIIIIIIIIIIIIIIIIIII	EF:f:0.005651
fffffffde3:GATGCGGGGGGGGGGGGGGG	16	chrX	639713	30	20M	*	0	0	CCCCCCCCCCCCCCCGCATC	IIIIIIIIIIIIIIIIIIII	EF:f:0.002585
ffffffff78:ATGCGGGGGGGGGGGGGGGG	16	chrX	639712	30	20M	*	0	0	CCCCCCCCCCCCCCCCGCAT	IIIIIIIIIIIIIIIIIIII	EF:f:0.008120

==> sacCer3.rep <==
fe095954fc:AGGGACCCCTCCCTAATGGG	16	chrM	84209	10	20M	*	0	0	CCCATTAGGGAGGGGTCCCT	IIIIIIIIIIIIIIIIIIII	EF:f:0.071630
fe095954fc:AGGGACCCCTCCCTAATGGG	0	chrM	85629	10	20M	*	0	0	AGGGACCCCTCCCTAATGGG	IIIIIIIIIIIIIIIIIIII	EF:f:0.071630
3f8256553f:GGGACCCCTCCCTAATGGGA	16	chrM	84208	10	20M	*	0	0	TCCCATTAGGGAGGGGTCCC	IIIIIIIIIIIIIIIIIIII	EF:f:0.305550
3f8256553f:GGGACCCCTCCCTAATGGGA	0	chrM	85630	10	20M	*	0	0	GGGACCCCTCCCTAATGGGA	IIIIIIIIIIIIIIIIIIII	EF:f:0.305550
cfe095954f:GGACCCCTCCCTAATGGGAG	16	chrM	84207	10	20M	*	0	0	CTCCCATTAGGGAGGGGTCC	IIIIIIIIIIIIIIIIIIII	EF:f:0.233040
cfe095954f:GGACCCCTCCCTAATGGGAG	0	chrM	85631	10	20M	*	0	0	GGACCCCTCCCTAATGGGAG	IIIIIIIIIIIIIIIIIIII	EF:f:0.233040
f3f8256553:GACCCCTCCCTAATGGGAGG	16	chrM	84206	10	20M	*	0	0	CCTCCCATTAGGGAGGGGTC	IIIIIIIIIIIIIIIIIIII	EF:f:0.195745
f3f8256553:GACCCCTCCCTAATGGGAGG	0	chrM	85632	10	20M	*	0	0	GACCCCTCCCTAATGGGAGG	IIIIIIIIIIIIIIIIIIII	EF:f:0.195745
aff0029a59:CTCCTTCTTAAAAAGGGGTT	16	chrM	85659	10	20M	*	0	0	AACCCCTTTTTAAGAAGGAG	IIIIIIIIIIIIIIIIIIII	EF:f:0.036516
30aa954353:GACCGAACCCCTTTTTAAGA	0	chrM	85654	10	20M	*	0	0	GACCGAACCCCTTTTTAAGA	IIIIIIIIIIIIIIIIIIII	EF:f:0.193518

