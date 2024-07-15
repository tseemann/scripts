#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fastq write_fasta assert);

my(@Options, $verbose, $fasta, $newid, $nons, $trunc, $nopoly);
setOptions();

my $nread=0;
my $nwrote=0;

while ( not eof() ) {
  my $l = read_fastq(\*ARGV);
  $nread++;

  # truncate first
  if ($trunc > 0) {
    $l->[1] = substr($l->[1], 0, $trunc);
    $l->[2] = substr($l->[2], 0, $trunc);
  }

  # then skip read if either left or right has N in it
  next if $nons and index($l->[1], 'N') >= $[ ;
  
  # don't allow monopolymers through ?
  next if $nopoly and ( $l->[1] =~ m/^(.)\1*$/i );
                       
  $nwrote++;

  if ($newid) {
    $l->[0] = "$nwrote/1";
  }

  if ($fasta) {
    write_fasta(\*STDOUT, $l);
  }
  else {
    write_fastq(\*STDOUT, $l);
  }
}

printf STDERR "Parsed %d reads. Discarded %d. Wrote %d (%.2f%%)\n", 
  $nread, 
  $nread-$nwrote, 
  $nwrote, ($nwrote*100/$nread);


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
    {OPT=>"n|nons!", VAR=>\$nons, DEFAULT=>0, DESC=>"Skip reads with Ns"},
    {OPT=>"p|nopoly!", VAR=>\$nopoly, DEFAULT=>0, DESC=>"Skip monopolymer reads"},
    {OPT=>"t|trunc=i", VAR=>\$trunc, DEFAULT=>0, DESC=>"Truncate reads to this length (0 = no truncation)"},
    {OPT=>"i|newid!", VAR=>\$newid, DEFAULT=>0, DESC=>"Generate new, compact IDs"},
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
  print "Usage: $0 [-fnpi] dirty.fq > clean.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
