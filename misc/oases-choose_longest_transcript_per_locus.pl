#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

# >Locus_1_Transcript_1/4_Confidence_0.250_Length_2051
# >Locus_1_Transcript_2/4_Confidence_0.667_Length_3861
# >Locus_1_Transcript_3/4_Confidence_0.750_Length_4131

my(@Options, $verbose, $minlen);
setOptions();

my $in  = Bio::SeqIO->new(-fh=>\*ARGV,   -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $locus = 0;
my $best;

while (my $seq = $in->next_seq) {
  print STDERR "\rProcessed: $." if $. % 9871 == 0;
  #                     $1               $2     $3                $4              $5
  $seq->id =~ m/^Locus_(\d+)_Transcript_(\d+)\/(\d+)_Confidence_([\d\.]+)_Length_(\d+)$/;
  if ($1 != $locus) {
    # print out best from last set 
    if ($locus and $best->length > $minlen) {
      print STDERR "$locus ", $best->id, " ", $best->length, "\n" if $verbose;
      $out->write_seq($best);
    }
    # start new set
    $best = $seq;
    $locus = $1;
  }
  elsif ($5 > $best->length) {
    # update best set
    $best = $seq;
  }
} 
print STDERR "\nDone.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"minlen=i",  VAR=>\$minlen, DEFAULT=>300, DESC=>"Minimum transcript length to keep"},
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
