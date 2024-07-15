#!/usr/bin/env perl
use strict;
use Data::Dumper;

#HWUSI-EAS-100R  0002    3       2       999     6334    0       1       ............................................................................    BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB    0

my(@Options, $verbose);
setOptions();

while (<ARGV>) {
  chomp;
  my @x = split m/\t/;
  next if substr($x[-1],0,1) eq '0'; # failed purity test
  print "\@", join(':',@x[0..7]),"\n";
  $x[8] =~ s/\./N/g; # replace . with N
  print $x[8],"\n+\n",$x[9],"\n";
}


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
  print "Usage: $0 [options] s_8_1_*_qseq.txt > s_8_1_sequence.fastq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
