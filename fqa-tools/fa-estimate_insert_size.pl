#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fasta);
use List::Util qw(sum min max);

# This script takes interleaved, fasta, pairs of reads
# and maps them to a reference fasta genome
# It then trims the reads to $trim bases (25 default)
# and looks for an EXACT MATCH in the reference genome \
# [OLD] using index()  
# [NEW] using a hash lookup of pre-computed $trim-mers
# if the left read hits the fwd strand, and right read the rev strand
# and they are separated less than $maxdep (500 default)
# the distance is recorded. Once we have $number (100 default) good
# distance samples, we print out the statistics for them.
# REMINDER:
# - only looks at fwd-> <-rev matches, not rev-> <-fwd matches!
# - only looks for EXACT matches
# - the $trim option to increase probability of exact match
# - it assumes interleaved + fasta format for reads
# - does not deal with multiple mapping in any way

BEGIN {
  $SIG{INT} = \&show_results;
}

my(@Options, $verbose, $maxsep, $number, $offset, $kmer);
setOptions();


# load the reference
my $rfile = shift @ARGV;
print STDERR "Loading reference: $rfile\n";
my $ref = Bio::SeqIO->new(-file=>$rfile, -format=>'fasta');
my $seq = $ref->next_seq;
print STDERR "Using sequence: ",$seq->display_id,"\n";
my $fwd = $seq->seq;
my $rev = $seq->revcom->seq;
my $len = $seq->length;
print STDERR "Sequence length: $len bp\n";

# make the hash lookup tables
my %fwdk;
my %revk;
print STDERR "Making $kmer-mer hashes from $rfile. Please wait...\n";
for my $i (0 .. $len-1-$kmer) {
  $fwdk{ substr($fwd,$i,$kmer) } ||= $i;
  $revk{ substr($rev,$i,$kmer) } ||= ($len - $i - 1);
}
printf STDERR "Forward $kmer-mer hash keys: %d\n", scalar keys %fwdk;
printf STDERR "Reverse $kmer-mer hash keys: %d\n", scalar keys %revk;

print STDERR "Processing reads: @ARGV\n";
my $p = 0;
my @dist;
my ($l, $r);

while (not eof () ) {         # please read "perldoc -f eof" !
  printf STDERR "\rProcessed %d pairs, %d matches... (CTRL-C for partial results)", ++$p, scalar(@dist);

  # read the mates
  my $L = read_fasta(\*ARGV);  
  my $R = read_fasta(\*ARGV);

  # choose the $offset..$kmer as our 'representative' for the left read
  $L->[1] = substr($L->[1], $offset, $kmer);

  # check if it has an exact match in the reference (don't continue otherwise)
  $l = $fwdk{$L->[1]} or next;

  # choose the right read rep
  $R->[1] = substr($R->[1], $offset, $kmer);

  # check if it is on the reverse strand
  $r = $revk{$R->[1]} or next;

  # ensure they are not wacky pairs
  next unless $r - $l <= $maxsep;

  # keep the distance in an array
  push @dist, ($r-$l);
  printf STDERR "\tmatch %d $l $r %d\n", scalar(@dist), $r-$l if $verbose;

  # stop if we have enough samples
  last if @dist >= $number;
}

show_results();
exit;

sub show_results {
  printf STDERR "\r\nDone - %d matches", 0+@dist;
  # print statistics
  exit -1 unless @dist > 0;
  my $mean = sum(@dist)/@dist;
  printf "mean = %d\n", $mean; 
  printf "sdev = %.1f\n", sqrt( sum( map { ($_ - $mean)**2 } @dist ) / @dist );
  my $median = $dist[int(@dist/2)];
  printf "median = %d\n", $median;
  printf "mad = %.1f\n", sum( map { abs($_ - $median) } @dist ) / @dist;
  printf "min = %d\n", min(@dist);
  printf "max = %d\n", max(@dist);
  exit;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"m|maxsep=i",  VAR=>\$maxsep, DEFAULT=>1000, DESC=>"Maximum separation to include in calculations"},
    {OPT=>"n|num=i",  VAR=>\$number, DEFAULT=>500, DESC=>"Number of valid read pairs to stop at"},
    {OPT=>"o|offset=i",  VAR=>\$offset, DEFAULT=>1, DESC=>"Use k-mer this many bases into read (0=true start)"},
    {OPT=>"k|kmer=i",  VAR=>\$kmer, DEFAULT=>25, DESC=>"Use this many bases from --offset as k-mer"},
  );

  (@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] reference.fasta interleaved-reads.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
