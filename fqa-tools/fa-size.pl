#!/usr/bin/env perl
use strict;
use warnings;

my $nseq = 0;
my $bp   = 0;

print STDERR "(reading from STDIN)\n" if not @ARGV;

while (<ARGV>) {
  if (m/>/) {
    $nseq++;
    next;
  }
  else {
    chomp;
    $bp += length $_;
  }
}

print "$bp letters in $nseq sequences\n";
