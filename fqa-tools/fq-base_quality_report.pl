#!/usr/bin/env perl
use strict;
#use List::Util qw(min max sum);
use List::MoreUtils qw(pairwise);
use Data::Dumper;

my(@Options, $verbose, $solexa, $maxnum);
setOptions();

# This script is designed to be fast!
# It does not validate the input .fq sequences
# It assumes all the reads are the same length
#
# Torsten Seemann, 6 May 2010

$a = $b; # hack to disable warnings about "main::a" used only once

if ($verbose) {
  print STDERR "Reading from: ", (@ARGV ? "@ARGV" : "(stdin)"), "\n";
}

my $N=0;  # number of reads
my $offset = $solexa ? 59 : 64; # ASCII offset for Q string
my %base;

while (<>) {    # ID
  my $nt = <>;   # sequence
  chomp $nt;
  my @nt = unpack '(A1)*', $nt;
#  print "nt=",join('.', @nt)."\n";
  scalar(<>);   # ID2
  my $qs = <>;  # quality string
  chomp $qs;
  my @qs = map { ord($_)-$offset } (unpack '(A1)*', $qs);
  pairwise { $base{$a}{$b}++ } @nt, @qs;
  print STDERR "\rProcessing: $N" if ($N++ % 59 == 0);
  last if $maxnum > 0 and $N >= $maxnum;
}
print STDERR "\r";

print STDERR "Processed $N reads.\n";

print Dumper(\%base);

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"s|solexa!", VAR=>\$solexa, DEFAULT=>0, DESC=>"Interpret as 'old' quality values (Solexa, Illumina < 1.3)"},
    {OPT=>"n|num=i", VAR=>\$maxnum, DEFAULT=>0, DESC=>"Stop after this many reads (0=do all)"},
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
