#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fastq write_fasta assert);

my(@Options, $verbose, $fasta, $newid);
setOptions();

assert(@ARGV >= 2, "please provide two .fq files as input");
assert(-r $ARGV[0], "can not open left reads: $ARGV[0]");
assert(-r $ARGV[1], "can not open right reads: $ARGV[1]");

open my $L, '<', $ARGV[0];
open my $R, '<', $ARGV[1];

if ($verbose) {
  print STDERR "Interleaving $ARGV[0] + $ARGV[1]\n";
  print STDERR "Output is FAST", ($fasta ? 'A' : 'Q'), ".\n";
  print STDERR "Please be patient...\n";
}

my $idnum=0;
my $nread=0;
my $writer = $fasta ? \&write_fasta : \&write_fastq ;

while (not eof $L) {
  my $l = read_fastq($L);
  my $r = read_fastq($R);
  $nread++;
  print STDERR "\rProcessed: $nread" if $nread % 1E5 == 0;


  if ($newid) {
    $idnum++;
    $l->[0] = "R$idnum/1";
    $r->[0] = "R$idnum/2";
  }
 
  $writer->(\*STDOUT, $l);
  $writer->(\*STDOUT, $r);
}


printf STDERR "\rParsed $nread read pairs..\n", $nread;


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"i|newid!", VAR=>\$newid, DEFAULT=>0, DESC=>"Generate new, compact IDs"},
    {OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
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
  print "Usage: $0 [--fasta] left.fq right.fq > interleaved.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
