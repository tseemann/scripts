#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $table, $meth, $chop, $clean);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

while (my $seq = $in->next_seq) {
  print STDERR "Translating: ", $seq->id, "\n" if $verbose;
  $seq = $seq->translate(-codontable_id=>$table);
  my $aa = $seq->seq;
  if ($meth) {
    substr $aa, 0, 1, 'M';
  }
  if ($chop) {
    $aa =~ s/\*$//;
  }
  if ($clean) {
    $aa =~ s/[^A-Z]/X/gi;
  }
  $seq->seq($aa);
  $out->write_seq($seq); 
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"table=i",  VAR=>\$table, DEFAULT=>11, DESC=>"Translation table"},
    {OPT=>"meth!",  VAR=>\$meth, DEFAULT=>0, DESC=>"Force Met-start codon "},
    {OPT=>"chop!",  VAR=>\$chop, DEFAULT=>0, DESC=>"Remove trailing stop codons"},
    {OPT=>"clean!",  VAR=>\$clean, DEFAULT=>0, DESC=>"Replace non-AA with X"},
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
  print "Usage: $0 [options] codingseq.fna > protein.faa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
