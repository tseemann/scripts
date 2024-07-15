#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/";
use Prokka qw(inform read_gff3 sort_gff3 write_gff3);

my(@Options, $verbose);
setOptions();

inform "Sorting GFF3 file(s): @ARGV";
my $gff = read_gff3(@ARGV);
sort_gff3($gff);
write_gff3(\*STDOUT, $gff);

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
  print "Usage: $0 [options] file.gff3 > sorted.gff3\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

