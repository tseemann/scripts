#!/usr/bin/env perl
use strict;
use List::Util qw(min max);

my(@Options, $verbose, $fasta, $reverse, $sort, $N);
setOptions();

# This script is designed to be fast!
# It does not validate the input .fq sequences
# Torsten Seemann, March 2009

my $kept=0;
my $total=0;
my %read;

print STDERR "Reading .fq ...\n";

while (my $id = <ARGV>) {
  substr($id,0,1,'>') if $fasta;  
  my $seq = uc(<ARGV>);
  scalar(<ARGV>);
  my $qual = <ARGV>;
  push @{$read{$seq}} , [ $id, $qual ];
  $total++;
}

print STDERR "Writing .fq ...\n";

for my $seq (keys %read) {
  my $list = $read{$seq};
  next unless @$list > $N;
  for my $r (@$list) {
    print $r->[0], $seq;
    print "+\n", $r->[1] unless $fasta;
    $kept++;
  } 
}

my $rare = $total - $kept;
printf STDERR "%d reads, %d duplicates (%.2f%%), %d $N-singletons (%.2f%%)\n",
  $total,
  $kept,
  $kept*100/$total, 
  $rare,
  $rare*100/$total;

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"N|num=i", VAR=>\$N, DEFAULT=>1, DESC=>"Up to this many dupes is a 'singleton'"},
    {OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
#    {OPT=>"r|reverse!", VAR=>\$reverse, DEFAULT=>0, DESC=>"Keep singletons, not duplicates"},
    {OPT=>"s|sort!", VAR=>\$sort, DEFAULT=>0, DESC=>"Sort output"}, );


  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] Reads.fq moreReads.fq ... > cleanReads.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__
; 59 -5 0.7597
< 60 -4 0.7153
= 61 -3 0.6661
> 62 -2 0.6131
? 63 -1 0.5573
@ 64 0 0.5000
A 65 1 0.4427
B 66 2 0.3869
C 67 3 0.3339
D 68 4 0.2847
E 69 5 0.2403
F 70 6 0.2008
G 71 7 0.1663
H 72 8 0.1368
I 73 9 0.1118
J 74 10 0.0909
K 75 11 0.0736
L 76 12 0.0594
M 77 13 0.0477
N 78 14 0.0383
O 79 15 0.0307
P 80 16 0.0245
Q 81 17 0.0196
R 82 18 0.0156
S 83 19 0.0124
T 84 20 0.0099
U 85 21 0.0079
V 86 22 0.0063
W 87 23 0.0050
X 88 24 0.0040
Y 89 25 0.0032
Z 90 26 0.0025
[ 91 27 0.0020
\ 92 28 0.0016
] 93 29 0.0013
^ 94 30 0.0010
_ 95 31 0.0008
` 96 32 0.0006
a 97 33 0.0005
b 98 34 0.0004
c 99 35 0.0003
d 100 36 0.0003
e 101 37 0.0002
f 102 38 0.0002
g 103 39 0.0001
h 104 40 0.0001 
