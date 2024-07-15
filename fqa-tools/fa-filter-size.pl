#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $minsize, $maxsize);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $nread=0;
my $nwrote=0;

while (my $seq = $in->next_seq) {
  $nread++;
  if ($seq->length >= $minsize and $seq->length <= $maxsize) {
    $out->write_seq($seq);
    $nwrote++;
  }
  print STDERR "\rWrote $nwrote/$nread" if $nread%1000==0;
} 

print STDERR "\rRead $nread sequences, wrote $nwrote, skipped ", $nread-$nwrote,".\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"minsize=i",  VAR=>\$minsize, DEFAULT=>0, DESC=>"Minimum sequence length"},
    {OPT=>"maxsize=i",  VAR=>\$maxsize, DEFAULT=>1E10, DESC=>"Maximum sequence length"},
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
  print "Usage: $0 [options] all.fasta > some.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
