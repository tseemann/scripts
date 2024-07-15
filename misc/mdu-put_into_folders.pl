#!/usr/bin/env perl
use warnings;
use strict;
use File::Spec;
use Data::Dumper;

my(@Options, $verbose, $nonmdu);
setOptions();

my %dir;

my $pattern = $nonmdu ? qr/^([^_]+)/ : qr/^(\d{4}-\d+)/;

for my $fastq (@ARGV) {
  my(undef,undef,$fname) = File::Spec->splitpath($fastq);
  print STDERR "Found: $fname\n";
  next unless $fname =~ $pattern;
  print STDERR "Matched: $1\n";
  push @{$dir{$1}}, $fastq;
}

print STDERR Dumper(\%dir);

for my $d (sort keys %dir) {
  print "mkdir -v \Q$d\E\n";
  for my $f (@{$dir{$d}}) {
    print "mv -v \Q$f\E \Q$d\E\n";
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"nonmdu!",  VAR=>\$nonmdu, DEFAULT=>0, DESC=>"Handle non-MDU IDs"},
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
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
