#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use List::Util qw(min max);

my(@Options, $verbose, $qval, $perline);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');

while (my $seq = $in->next_seq) {
  print ">", $seq->display_id, " ", $seq->desc, "\n";
  my $L = $seq->length;
  while ($L > 0) {
    my $N = min($L, $perline);
    my $string = "$qval "x$N;
    chop $string;
    print $string,"\n";
    $L = $L - $perline;
  }
} 

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"q|quality=i",  VAR=>\$qval, DEFAULT=>40, DESC=>"Default quality value"},
    {OPT=>"n|perline=i",  VAR=>\$perline, DEFAULT=>20, DESC=>"Number of quality integers per line"},
  );

  (!@ARGV) && (usage());

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
