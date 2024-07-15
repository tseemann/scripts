#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $qualdesc);
setOptions();

# TODO - use the read_fastq library

while (<ARGV>) {
  chomp;
  print ">",substr($_,1);
  my $seq = <ARGV>;
  $_ = <ARGV>;
  my $qual = <ARGV>;
  print ($qualdesc ? " $qual" : "\n");
  print $seq;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"qd|qualdesc!", VAR=>\$qualdesc, DEFAULT=>0, DESC=>"Put quality string as FASTA description"},
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
  print "Usage: $0 reads.fq [more.fq ...] > combined.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
