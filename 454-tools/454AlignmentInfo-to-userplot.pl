#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(max);

#Position              Reference             Consensus             Quality Score         Unique Depth          Align Depth           Total Depth           Signal                StdDeviation
#>gi|15644634|ref|NC_000915.1|               154
#154                   G                     G                     11                    1                     1                     1                     1.02                  1.02
#155                   T                     T                     11                    1                     1                     1                     0.64                  0.64
#156                   G                     -                     1                     1                     1                     1                     0.00                  0.00

my(@Options, $verbose);
setOptions();

my $header  = scalar <ARGV>;
die "Does not look like 454AlignmentInfo.tsv file" unless $header =~ m/^Position/;

my %freq;

while (<ARGV>) {
  chomp;
  my @x = split /\t/;
  next unless @x >= 8;
  my $pos = $x[0];
  next unless $pos =~ m/\d+/;
  $freq{$pos} = $x[6];
  print STDERR "\rReading: pos" if $pos%789==0;
}

my $L = max(keys %freq);
print STDERR "\rWriting $L lines...\n";
for my $i (1 .. $L) {
  my $c = $freq{$i} || 0;
  print "$c\n";
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
  print "Usage: $0 [options] 454AlignmentInfo.tsv > 454.userplot\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
