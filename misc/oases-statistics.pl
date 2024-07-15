#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

# >Locus_1_Transcript_1/4_Confidence_0.250_Length_2051
# >Locus_1_Transcript_2/4_Confidence_0.667_Length_3861
# >Locus_1_Transcript_3/4_Confidence_0.750_Length_4131

my(@Options, $verbose, $minlen);
setOptions();

my $count=0;
my $uniq=0;
my %locus;
my @len;
my @conf;

while (<>) {
  #                       1                2      3                 4               5
  next unless m/^>Locus_(\d+)_Transcript_(\d+)\/(\d+)_Confidence_([\d\.]+)_Length_(\d+)$/;
  my $i = $1;
  my $c = $4;
  my $L = $5;
  $uniq++ if $3 == 1;
  $count++;
  $locus{$i}++;
  $len[int($L/1000)]++;
  $conf[int($c*10)]++;
}

printf "Transcripts\t%d\n", $count;
printf "Unique Locii\t%d\n", $uniq;
printf "Locii\t%d\n", scalar keys %locus;

for my $i (0 .. $#len) {
  if ($len[$i]) {
    printf "Length:%d-%d\t%d\n", 1000*$i, 1000*($i+1)-1, $len[$i];
  }
}

for my $i (0 .. $#conf) {
  if ($conf[$i]) {
    printf "Conf:%d-%d%%\t%d\n", $i*10, 10*($i+1)-1, $conf[$i];
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"minlen=i",  VAR=>\$minlen, DEFAULT=>200, DESC=>"Minimum transcript length to keep"},
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
  print "Usage: $0 [options] oases_transcripts.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
