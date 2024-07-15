#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $reference, $coverage);
setOptions();

my $refsize = 1E20;  
if ($reference =~ m/^(\d+)$/) {
  $refsize = $1;
}
elsif (-r $reference) { 
  $refsize = -s $reference; # file size is ~ genome size
}
print STDERR "Target is ${coverage}x of $refsize bp genome.\n";
my $target = $refsize * $coverage;
my $nseq=0;
my $sum=0;

while (<>) {
  print;
  my $seq = <>;
  $sum += length($seq)-1;
  $nseq++;
  print $seq;
  print scalar(<>);
  print scalar(<>);
  if ($. % 4000 == 0) {  # no need to check after every read
    print STDERR "### nseq=$nseq sumlen=$sumlen\n" if $verbose;
    last if $sum >= $target;
  }
}
print STDERR "Wrote $nseq reads totalling $sum bp.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"reference=s", VAR=>\$reference, DEFAULT=>'', DESC=>"Reference genome (text) or size in bp (integer)"},
    {OPT=>"coverage=f", VAR=>\$coverage, DEFAULT=>25, DESC=>"Desired coverage"},
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
  print "Usage: $0 < reads.fq > subset.fq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
