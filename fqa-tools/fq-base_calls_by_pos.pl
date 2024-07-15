#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq int_qual);
use List::Util qw(sum);

my(@Options, $verbose, $update, $noheader, $nopc, $showfreq);
setOptions();

my %freq;
my $count=0;

while ( not eof() ) {
  my $read = read_fastq(\*ARGV);

#  SLOWER?
#  my $pos=0;
#  for my $base (split m//, $read->[1]) {
#    $freq{++$pos}{$base}++;
#  }

  # FASTER
  my $L = length($read->[1]);
  $freq{$_}{substr($read->[1],$_,1)}++ for (0 .. $L-1);

  $count++;
  print STDERR "\rProcessing: $count" if ($count % $update == 0);
}
print STDERR "\n";

print STDERR Dumper(\%freq) if $verbose;

my @BASES = qw(A T G C N);
print join("\t", 'POS', 'TOTAL', @BASES),"\n" unless $noheader;

for my $pos (sort { $a <=> $b } keys %freq) {
  print "$pos";
  my $total = sum( values %{$freq{$pos}} );
  print "\t$total";
  for my $base (@BASES) {
    my $f = $freq{$pos}{$base} || 0;
    print "\t$f" if $showfreq;
    printf("\t%6.2f%%", $f*100/$total) unless $nopc;
  }
  print "\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"u|update=i", VAR=>\$update, DEFAULT=>10000, DESC=>"Update frequency"},
    {OPT=>"noheader!", VAR=>\$noheader, DEFAULT=>0, DESC=>"Don't print header"},
    {OPT=>"nopc!", VAR=>\$nopc, DEFAULT=>0, DESC=>"Don't print percents"},
    {OPT=>"freq!", VAR=>\$showfreq, DEFAULT=>0, DESC=>"Print raw frequencies"},
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
  print "Usage: $0 [<] reads.fq ...\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
