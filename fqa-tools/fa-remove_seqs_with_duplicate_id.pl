#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose);
setOptions();

my $nread=0;
my $nwrote=0;
my %seen;
my $keep=1;

while (<ARGV>) {
  if (m/^>(\S+)/) {
    $keep = not $seen{$1};
    $seen{$1}++;
    $nread += 1;
    $nwrote += $keep;
    print STDERR "\rProcessed: $nread" if $verbose and $nread%10000==0;
    print STDERR "\nRemoving duplicate: $_" if !$keep and $verbose;
  }
  print if $keep;  
}

print STDERR "\nRead $nread, ",$nread-$nwrote, " duplicates, wrote $nwrote.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose info"},
#    {OPT=>"u|unique!",  VAR=>\$unique, DEFAULT=>0, DESC=>"Don't even output one representative for duplicates"},
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
  print "Usage: $0 [options] seq.fasta [ more.fasta ... ] > unique.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
