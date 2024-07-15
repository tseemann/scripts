#!/usr/bin/env perl
use strict;
use warnings;

my $ss = shift @ARGV;
$ss ||= 'SampleSheet.csv';
open my $SS, '<', $ss;
LINE:
while (<$SS>) {
  my @x = split m/,/;
  my $lane = $x[1];
  next unless $lane =~ m/^[12345678]$/;
  my $name = $x[2];
  $name =~ s{/}{_}g;
  my $barcode = $x[4];
  next unless $barcode =~ m/^[AGTC]{6}$/;
  for my $dir ('1', '2') {
    my $sf = "s_${lane}_${dir}_${barcode}.txt";
    next LINE unless -r $sf;
    print "mkdir -p $name/analysis\n";
    print "mv -v s_${lane}_${dir}_${barcode}.txt $name/\n";
  }
  print "fq-pre-process-pe -a $name/*.txt > $name/analysis/s_paired.fa\n";
}

__DATA__

FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator
635F4AAXX,1,Sal_1685,none,GCCAAT,264,N,NA,Scott_Coutts
635F4AAXX,1,EC3134,none,CAGATC,265,N,NA,Scott_Coutts
635F4AAXX,1,EC3204,none,ACTTGA,266,N,NA,Scott_Coutts
635F4AAXX,1,CS1,none,GATCAG,267,N,NA,Scott_Coutts
635F4AAXX,2,Bp_06872-2,none,GATCAG,272,N,NA,Scott_Coutts
635F4AAXX,2,Bp_WT-2,none,TAGCTT,273,N,NA,Scott_Coutts
635F4AAXX,3,Ab_6009-1,none,TAGCTT,231,N,NA,Scott_Coutts
635F4AAXX,3,Ab_6009-2,none,GGCTAC,232,N,NA,Scott_Coutts
635F4AAXX,3,Ab_1814,none,CTTGTA,233,N,NA,Scott_Coutts
635F4AAXX,3,Ab_1815,none,ATCACG,234,N,NA,Scott_Coutts
635F4AAXX,3,Ab_1811,none,CGATGT,235,N,NA,Scott_Coutts
635F4AAXX,3,Abces,none,TTAGGC,202,N,NA,Scott_Coutts
635F4AAXX,3,Hj5,none,TGACCA,203,N,NA,Scott_Coutts
635F4AAXX,4,PhiX174,PhiX174,TTAGGC,PhiX174,Y,NA,Scott_Coutts
635F4AAXX,5,Mu_000945,none,CGATGT,260,N,NA,Scott_Coutts
635F4AAXX,5,Mu_991845,none,TTAGGC,261,N,NA,Scott_Coutts
635F4AAXX,5,Mu_001506,none,TGACCA,262,N,NA,Scott_Coutts
635F4AAXX,5,Mu_980535,none,ACAGTG,263,N,NA,Scott_Coutts
635F4AAXX,6,Sa_Eur,none,ACAGTG,268,N,NA,Scott_Coutts
635F4AAXX,6,Sa_Taiwan,none,GCCAAT,269,N,NA,Scott_Coutts
635F4AAXX,6,Mm_MON17,none,CAGATC,270,N,NA,Scott_Coutts
635F4AAXX,6,Mm_MON4,none,ACTTGA,271,N,NA,Scott_Coutts
635F4AAXX,7,FA1090_PilE+AV,none,ATCACG,223,N,NA,Scott_Coutts
635F4AAXX,7,FA1090_PilE-AV,none,CGATGT,224,N,NA,Scott_Coutts
635F4AAXX,7,FA1090_Maf+AV,none,TTAGGC,226,N,NA,Scott_Coutts
635F4AAXX,8,FA1090_Opa+AV,none,TGACCA,227,N,NA,Scott_Coutts
635F4AAXX,8,FA1090_Opa-AV,none,ACAGTG,228,N,NA,Scott_Coutts
635F4AAXX,8,FAM18_PilE2_AV,none,GCCAAT,229,N,NA,Scott_Coutts
635F4AAXX,8,NMB_PilE2_AV,none,CAGATC,230,N,NA,Scott_Coutts
