#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fasta);

my(@Options, $verbose);
setOptions();

#my $writer_func = $fasta ? \&write_fasta : \&write_fastq;

my $npair=0;

while (not eof () ) {         # please read "perldoc -f eof" !
  my $L = read_fasta(\*ARGV);
  my $R = read_fasta(\*ARGV);
  print ">", $L->[0], "\n", $L->[1], ":", $R->[1], "\n";
  print STDERR "\rProcessing pairs: $npair" if (++$npair % 9871 == 0);
}
print STDERR "\rWrote $npair pairs of reads in SSAKE format.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
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
  print "Usage: $0 [options] reads.fq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
