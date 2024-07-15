#!/usr/bin/perl -w
use strict;

my(@Options, $inverse);
setOptions();

#0 readname contigname strand 
#3 contigstart contigend 
#5 readstart readend readlength 
#8 score editstring

my $count=0;

while (my $line = <ARGV>) {
  chomp $line;
  my @f = split m/\t/, $line;
  next unless @f == 10;
  my $exact = ($f[5] == 1 and $f[6] == $f[7] and $f[9] eq $f[7]);
  print $line,"\n" if $exact && !$inverse or !$exact && $inverse;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v!",  VAR=>\$inverse, DEFAULT=>0, DESC=>"Return only NON-EXACT hits instead"},
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
