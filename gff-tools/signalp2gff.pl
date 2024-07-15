#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use List::Util qw(min max);

use FindBin;
use lib "$FindBin::Bin/";
use Prokka qw(inform make_id);

## SignalP-HMM gram- predictions
## name      !  Cmax  pos ?  Sprob ?
#L550/2/13786/12461  Q  0.000   0 N  0.000 N  
#L550/2/14691/13786  Q  0.257  25 N  0.288 N  
#L550/2/17743/16394  Q  0.000   0 N  0.000 N  
#L550/2/18793/17843  S  0.474  25 N  0.643 Y  
#L550/2/20101/18800  Q  0.000  28 N  0.000 N  
#L550/2/20091/21191  Q  0.000  17 N  0.000 N  
#L550/2/21298/22404  S  0.509  29 Y  0.775 Y  
#L550/2/23525/22608  Q  0.000   0 N  0.000 N  
#L550/2/26735/26478  S  0.949  23 Y  0.978 Y  

my(@Options, $verbose, $gff_version, $source);
setOptions();

my $count=0;

my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gff_version);
my $f = Bio::SeqFeature::Generic->new;
$f->gff_format($gff_factory);


while (<ARGV>) {
  chomp;
  my @x = split m/\t/;
  next unless @x == 9 and $x[1] =~ m/^\d+$/ and $x[5] =~ m/^[ATCG]{3}$/i;

  # strip spaces from either end (the coordinates seem to have post-spaces)
  map { s/^\s+// } @x;
  map { s/\s+$// } @x;

  my $ID = make_id('trna', ++$count);
  
  inform "[$count] Found $ID $x[4]($x[5]) in $x[0]\n" if $verbose;

  my $f = Bio::SeqFeature::Generic->new( 
    -seq_id     => $x[0],
    -source_tag => 'SignalP-3.0',
    -primary    => 'signal_peptide',
    -start      => min($x[2], $x[3]),
    -end        => max($x[2], $x[3]),
    -strand     => ($x[2] <= $x[3] ? +1 : -1),
    -tag        => {
      'ID' => $ID,
      'Parent' => $ID,
      'product' => "signal peptide",
    }
  );

  print "##gff-version $gff_version\n" if $count==1;

  print $f->gff_string($gff_factory), "\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"gff=i", VAR=>\$gff_version, DEFAULT=>3, DESC=>"GFF Version: 1, 2, 3"},
    {OPT=>"source=s", VAR=>\$source, DEFAULT=>'tRNAScan', DESC=>"tRNA prediction software"},
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
  print "Usage: $0 [options] < orfs.fsa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
