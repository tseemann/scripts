#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $num);

# try and emulate tradtional Unix 'head' syntax: % head -10 file.txt
BEGIN {
  if (@ARGV and $ARGV[0] =~ m/^-(\d+)$/ and $1 >= 0) {
    $num = $1;
    shift @ARGV;
  }
}

setOptions();

my $in  = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
while ($num-- > 0) {
  my $seq = $in->next_seq or last;
  $out->write_seq($seq);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"verbose info"},
    {OPT=>"n|num=i",  VAR=>\$num, DEFAULT=>10, DESC=>'Write the first N sequences out'},
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
