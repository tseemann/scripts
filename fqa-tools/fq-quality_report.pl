#!/usr/bin/env perl
use strict;
use List::Util qw(min max sum);
use List::MoreUtils qw(first_index last_index);
use Data::Dumper;

my(@Options, $verbose, $solexa, $horizontal);
setOptions();

# This script is designed to be fast!
# It does not validate the input .fq sequences
# It assumes all the reads are the same length
#
# Torsten Seemann, 30 Mar 2009

if ($verbose) {
  print STDERR "Reading from: ", (@ARGV ? "@ARGV" : "(stdin)"), "\n";
}

my @qual;
my $N=0;  # number of reads
my $offset = $solexa ? 59 : 64; # ASCII offset for Q string

while (<>) {    # ID
  scalar(<>);   # sequence
  scalar(<>);   # ID2
  my $qs = <>;  # quality string
#  print $qs;
  chomp $qs;
  my $L = length($qs);
  for my $i (0 .. $L-1) {
    my $ql = substr($qs,$i,1);
    my $q = ord($ql) - $offset;  # http://en.wikipedia.org/wiki/FASTQ_format
    $qual[$i][$q]++;
  }
  $N++;
  print STDERR "\rProcessing: $N" if ($N % 98765 == 0);
}
print STDERR "\r";

print STDERR "Processed $N reads.\n";
my $L = scalar(@qual);
print STDERR "Assuming read length is $L.\n";

if ($horizontal) {
  print "POSITION ", join(" ", map { sprintf("%02d",$_) } (1..$L)),"\n";
  print "QUALITY  ";
  for my $i (0 .. $L-1) {
# Illumina 1.0
#    my $wm = sum( map { ($qual[$i][$_+5] || 0) * $_ } (-5..40) ) / $N;
# Illumina 1.3+ http://en.wikipedia.org/wiki/FASTQ_format
    my $wm = sum( map { ($qual[$i][$_] || 0) * $_ } (0..62) ) / $N;
    printf "%02d ",$wm;
  }
  print "\n";
}
else {
  print "POS\tMIN\tAVG\tMAX\tSYM\n";
  my $twm=0;
  my $min=100;
  my $max=-100;
  for my $i (0 .. $L-1) {
    my $minq = first_index { $_ } @{$qual[$i]};
    my $maxq = last_index  { $_ } @{$qual[$i]};
    my $wm = sum( map { ($qual[$i][$_] || 0) * $_ } (0..62) ) / $N;
    printf "%d\t%d\t%4.1f\t%d\t%c\n", $i+1, $minq, $wm, $maxq, int(64+$wm+0.5);
    $twm+=$wm;
    $min = min($min,$minq);
    $max = max($max,$maxq);
  }
  printf "MEAN\t%d\t%4.1f\t%d\t\#\n", $min, $twm/$L, $max;
}

#print Dumper(\@qual);

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"s|solexa!", VAR=>\$solexa, DEFAULT=>0, DESC=>"Interpret as 'old' quality values (Solexa, Illumina < 1.3)"},
    {OPT=>"z|horizontal!", VAR=>\$horizontal, DEFAULT=>0, DESC=>"Horizontal report"},
  );

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] < solexa.fastq\n";
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
