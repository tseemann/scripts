#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

my(@Options, $verbose);
setOptions();

my $in  = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $old = 0;
my $new = 0;

while (my $seq = $in->next_seq) {
  $old++;
  my $s = $seq->id.' '.$seq->desc;
  my @d = split m/\001/, $s;        # CTRL-A is \001
  for my $d (@d) {
    print ">$d\n", $seq->seq, "\n";
    $new++;
  }
} 

print STDERR "Read $old, wrote $new sequences after CTRL-A separating.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
