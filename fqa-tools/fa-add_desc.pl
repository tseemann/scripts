#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $desc, $append);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

while (my $seq = $in->next_seq) {
  print STDERR "\rProcessing: $seq->id";
  my $new = $append ? $seq->desc.$desc : $desc;
  $seq->desc($new);
  $out->write_seq($seq);
} 
print STDERR "\rDone.                                     \n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"d|desc=s",  VAR=>\$desc, DEFAULT=>'', DESC=>"Set all sequences to this description"},
    {OPT=>"a|append!",  VAR=>\$append, DEFAULT=>0, DESC=>"Append to existing desc, don't replace"},
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
  print "Usage: $0 [options] old.fasta > new_with_descriptions.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
