#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq);
use List::Util qw(sum);

use constant THROBBER => 9751;

$SIG{INT} = \&show_result;

my(@Options, $verbose);
setOptions();

my $count=0;
my %freq;

while ( not eof() ) {
  $count++;
  print STDERR "\rProcessing: $count (CTRL-C for partial results)" if ($count % THROBBER == 0);
  my $read = read_fastq(\*ARGV);
  my $geecee = ($read->[1] =~ tr/GC/GC/);
  $freq{'gc'} += $geecee;
  $freq{'bp'} += length($read->[1]);
  next unless $read->[1] =~ m/GGC(.)G/;
  $freq{$1}++;
}
show_result();


sub show_result {
  print STDERR "\r", ' 'x50, "\r";
  print "Processed $count reads.\n";
  printf "GC = %.2f%%\n", $freq{gc}*100/$freq{bp};
  print Dumper(\%freq);
  exit;
}  

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
#    {OPT=>"u|update=i", VAR=>\$update, DEFAULT=>9751, DESC=>"Update frequency"},
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
  print "Usage: $0 [<] reads.fq ...\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
