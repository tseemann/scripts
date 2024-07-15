#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fastq assert write_fasta);

my(@Options, $verbose, $fasta, $left, $right);
setOptions();

my $fname = "@ARGV" || '(stdin)';
my $writer_func = $fasta ? \&write_fasta : \&write_fastq;

while (not eof () ) {         # please read "perldoc -f eof" !
  my $read = read_fastq(\*ARGV);
  my $L = length($read->[1]);
  for my $i (1, 2) {
    $read->[$i] = substr($read->[$i], $left, $L-$right-$left);
  }
  $writer_func->(\*STDOUT, $read);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
    {OPT=>"l|left=i", VAR=>\$left, DEFAULT=>0, DESC=>"Trim first X bases"},
    {OPT=>"r|right=i", VAR=>\$right, DEFAULT=>0, DESC=>"Trim last X bases"},
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
  print "Usage: $0 [options] long.fq > shorter.fq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
