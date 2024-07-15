#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Data::Dumper;

my(@Options, $verbose);
setOptions();

#----------------------------------------------------------------------

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');

my %freq;

while (my $seq = $in->next_seq) {
  for my $letter (qw(A T C G)) {
    print STDERR $seq->id,"\t$letter\n";
    count_runs( uc $seq->seq, $letter );
  }
}
print Dumper(\%freq); 

#----------------------------------------------------------------------

sub count_runs {
  my($dna, $letter) = @_;
  # pos($s) is the position (base-0) of the first non-char AFTER run (match)
  # but remember DNA is (base-1) ... sigh
  while ($dna =~ m/($letter+)/gi) {
    my $L = length $1;
    $freq{$L}{$letter}++;
    pos($dna) += $L;
  }  
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
#    {OPT=>"chars=s",  VAR=>\$chars, DEFAULT=>'N', DESC=>"Characters to consider in runs (case-insenstive)"},
#    {OPT=>"flank=i",  VAR=>\$flank, DEFAULT=>16, DESC=>"Flanking length to show"},
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
  print "Usage: $0 [options] long_scaffolds.fasta > shorter.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
