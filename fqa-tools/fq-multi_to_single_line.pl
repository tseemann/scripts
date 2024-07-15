#!/usr/bin/env perl
use strict;

my(@Options, $verbose, $qualdesc);
setOptions();

my($id,$seq,$qual,$wantq);

while (<ARGV>) {
  chomp;
  if (m/^\@/) {
    print "$id\n$seq\n+\n$qual\n" if $wantq;
    $id = $_;
    $seq = $qual = $wantq = undef;
  }
  elsif (m/^\+/) {
    $wantq=1;
    next;
  }
  elsif ($wantq) {
    $qual .= $_;
  }
  else {
    $seq .= $_;
  }
}

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
  print "Usage: $0 wrapped.fq [more_wrapped.fq ...] > singlelines.fq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
