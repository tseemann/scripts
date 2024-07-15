#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use List::Util qw(min max);

use FindBin;
use lib "$FindBin::Bin/";
use Prokka qw(inform make_id);

my(@Options, $verbose, $gff_version, $source);
setOptions();

#Tandem Repeats Finder Program writen by: 
#Version 4.00
#Sequence: pMUM003/1 M. marinum DL240490 plasmid ; JKD2945S010002H03R DL.5.563 
#Parameters: 2 7 7 80 10 50 500
#1274 1335 26 2.3 26 81 10 63 50 9 3 37 1.52 TATATAATACAATAAGAATATAATAT TATATAACTCCAATAAGATTATAAATATTATATAATACAATAG AATATACTATTATCATAAT
#28350 28375 11 2.4 11 100 0 52 76 0 7 15 0.99 AAAATGAATAA AAAATGAATAAAAAATGAATAAAAAA
#37875 37956 42 2.0 42 82 0 101 36 15 23 24 1.94 TAGATGCTAAAGCGTCAGGAATAACATCAATAGCAATAGGAG TAGATTCTAAGGCTTCAGGACTT ACATCAATGGCAATAGGAGTAGATGCTAAAGCGTCAGGAATAACTTCAATAGCAATAGG
#38804 38861 21 2.6 21 83 10 71 67 1 13 17 1.32 ATAAAAAAATAGAAAAAATTG ATAAAAGAAATAGACAAAAAGTTAGATGAAAAAATAGAAAAAATT GATAAAAAAATAG

my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gff_version);

my $seqid='';
my $version='';
my $count=0;

while (<ARGV>) {
  chomp;
  if (m/^Sequence:\s+(\S+)/) {
    $seqid = $1;
    inform "Examining: $seqid\n" if $verbose;
  }
  elsif (m/^Version\s+(\S+)/) {
    $version = $1;
    inform "Detected TRF version '$version'\n" if $verbose;
  }
  elsif (m/^\d+\s/) {
    my @x = split m/\s+/;
    next unless @x == 15;

    if (length($x[13]) > 30) {
      $x[13] = substr($x[13],0,14).'..'.substr($x[13],-14,14);
    }

    my $ID = make_id('trf', ++$count);

    inform "[$count] Found $ID $x[4]($x[5]) in $seqid\n" if $verbose;

    my $f = Bio::SeqFeature::Generic->new( 
      -seq_id     => $seqid,
      -source_tag => $source,
      -primary    => 'tandem_repeat',
      -start      => min($x[0], $x[1]),
      -end        => max($x[0], $x[1]),
      -strand     => ($x[0] <= $x[1] ? +1 : -1),
#      -score      => $x[8],
      -tag        => {
	'ID' => $ID,
	'product' => "tandem repeat: $x[3] units x $x[2] bp ($x[13])",
      }
    );

    print "##gff-version $gff_version\n" if $count==1;

    print $f->gff_string($gff_factory), "\n";
  }
}
inform "Found $count features.\n" if $verbose;

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"gff=i", VAR=>\$gff_version, DEFAULT=>3, DESC=>"GFF Version: 1, 2, 3"},
    {OPT=>"source=s", VAR=>\$source, DEFAULT=>'trf', DESC=>"Prediction software"},
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
