#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose);
setOptions();

my $in  = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');

my $gc = 0;
my $tot = 0;

while (my $seq = $in->next_seq) {
  $tot += $seq->length;
  my $s = $seq->seq;
  $s =~ s/[^gGcC]//g;
  $gc += length $s;
}

printf "GC\t$gc\t$tot\t%.2f\n", $gc*100/$tot;

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"verbose info"},
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
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
