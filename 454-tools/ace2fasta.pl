#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my(@Options, $debug, $reads, $contigs);
setOptions();

die "please specify one of --reads or --contigs only" unless $reads xor $contigs;
my $tag = $reads ? 'RD' : 'CO';

my $inside=0;

while (<ARGV>) {
  chomp;
  if (m/^$tag\s+(.*)$/) {
    print ">$1\n";
    $inside=1;
  }
  elsif ($inside and m/^$/) {
    $inside=0;
  }
  elsif ($inside) {
    s/\*//g;
    s/[^AGTC]/n/gi;
    print "$_\n";
  }
  else {
    next;
  }
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"contigs!",  VAR=>\$contigs, DEFAULT=>0, DESC=>"Extract 'CO' seqs (contigs)"},
    {OPT=>"reads!",  VAR=>\$reads, DEFAULT=>0, DESC=>"Extract 'RD' seqs (reads)"},
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
