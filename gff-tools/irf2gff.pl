#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use List::Util qw(min max);

use FindBin;
use lib "$FindBin::Bin/";
use Prokka qw(inform make_id);

my(@Options, $verbose, $gff_version, $source, $minid, $minlen, $maxsep);
setOptions();

#Inverted Repeats Finder Program writen by:
#Version 3.05
#Sequence: S_moni/1
#Parameters: 2 3 5 80 10 40 100000 500000 -t7 600000 -d -h 

#911 993 83 1208 1289 82 214 73.2558 8.1395 43 89.0909 10.9091 95.2381 4.7619 0.0000 2201 2200
#TTATTTAAATTA AATATTGA ATATATTATTATCAATATTAAGCTTATATCTATTTTTTTATTAAAACCTATAAAAATAACTAT
# ATAGTTTTACTATCATTTTTTAAAAAAATATATG ATATGTACATATATAT CAAATAAAATATTAAATATATAACTCCAATAA

#28104 28183 80 28275 28354 80 91 83.7500 0.0000 95 54.3750 45.6250 53.7313 46.2687 0.0000 56458 56458 TTTCT
#TAACGCCAAGTTAGCAGACTACGCACTGCGTGCTTCGCTGGCGTAACTTGGCTAAAATTTTCCTTTCCCTTCGGG CCCGAATGAAGAGCTTATTTTTAGAAAAGTT
#ACGTCAGCGAAGCACGCAGTGCGTAGTCTGCTAACTTTTCGTGAAGAAA

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
    next unless @x == 19;

    next unless $x[7] >= $minid;
    next unless $x[2] >= $minlen or $x[5] >= $minlen;
    next unless $x[6] <= $maxsep;
    
    # stupid joiner sequence from TIGR for S_moni genome
    next if $x[16] eq 'CACACACTTAATTAATTAAGTGTGTG';
    next if $x[17] eq 'CACACACTTAATTAATTAAGTGTGTG';
    
    $x[7] = int($x[7]); # %id 
    $x[8] = int($x[8]); # %indel (%gap)

    my $ID = make_id('irf', ++$count);

    inform "[$count] Found $ID in $seqid\n" if $verbose;

    print "##gff-version $gff_version\n" if $count==1;

    for my $right (0, 1) 
    {
      my $dir = $right ? 'right' : 'left';
      my $dirl = $right ? 'L' : 'R';
      my $len = $right ? $x[5] : $x[2];
      my($begin,$end) = $right ? ($x[3],$x[4]) : ($x[0],$x[1]); # is this correct?
          
      my $f = Bio::SeqFeature::Generic->new( 
	-seq_id     => $seqid,
	-source_tag => $source,
	-primary    => 'inverted_repeat',
	-start      => $begin,
	-end        => $end,
	-strand     => $right ? -1 : +1,
#	-score      => $x[9],
	-tag        => {
	  'ID' => "$ID.$dirl",
	  'product' => "inverted repeat ($dir): ${len}bp $x[7]%id $x[8]%gap",
	}
      );

      print $f->gff_string($gff_factory), "\n";
    }
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
    {OPT=>"minid=f", VAR=>\$minid, DEFAULT=>70, DESC=>"Minimum %identity between repeat units"},
    {OPT=>"minlen=i", VAR=>\$minlen, DEFAULT=>10, DESC=>"Minimum repeat unit length"},
    {OPT=>"maxsep=i", VAR=>\$maxsep, DEFAULT=>10_000, DESC=>"Maximum separation between repeat units"},
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
