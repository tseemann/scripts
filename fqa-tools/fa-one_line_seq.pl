#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');

while (my $seq = $in->next_seq) {
  print STDERR "\rProcessed: $." if $. % 9871 == 0;
  print ">", $seq->display_id;
  print " ",$seq->desc if $seq->desc;
  print "\n";
  print $seq->seq, "\n";
} 
print STDERR "\nDone.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
#    {OPT=>"w|width=i",  VAR=>\$width, DEFAULT=>0, DESC=>"Preferred width: 0=unlimited"},
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
  print "Usage: $0 [options] multi_line.fasta [ multi_line2.fasta ... ] > single_line.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
