#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use List::Util qw(max);
use List::MoreUtils qw(any);

my(@Options, $verbose, $mink, $maxk);
setOptions();
$|=1;

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');

while (my $seq = $in->next_seq) {
  my $dna = uc($seq->seq);
  print $seq->display_id,"\t",length($dna);
  $maxk = length($dna);
  my $k = 1;
  my $dup;
  while ($k++ <= $maxk) {
    printf "%10d", $k;
    $dup = has_duplicate_kmer($dna, $k);
    last if !defined $dup;
    print STDERR $seq->display_id,"\t$k\t$dup\n" if $verbose;
    print "\b"x10;
  }
  print "\t$k\t$dup\n";
}


sub has_duplicate_kmer {
  my($dna, $k) = @_;
  my $L = length($dna);
  my %num_of;
  for my $i (0 .. $L-$k) {
    my $kmer = substr($dna, $i, $k);
    return $kmer if ++$num_of{$kmer} > 1;
  }
  return;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
#    {OPT=>"mink", VAR=>\$mink, DEFAULT=>1, DESC=>"Minimum k to test (0 = 1)"},
    {OPT=>"maxk", VAR=>\$maxk, DEFAULT=>0, DESC=>"Maximum k to test (0 = sequence length)"},
  );

#  (@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [--fasta] left.fq right.fq > interleaved.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
