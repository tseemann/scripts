#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;

use FindBin;
use lib "$FindBin::Bin/";
use Prokka qw(inform make_id);

my(@Options, $verbose, $gff_version, $source);
setOptions();

my $count=0;
my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gff_version);
my $seqid;
while (<ARGV>) {
  chomp;
  if($_ =~ m/^>(\S+)/){
	  $seqid = $1;
	  next;
  }
  my @x = split m/\s+/;
  next unless @x == 5 and $x[0] =~ m/^\d+$/ and $x[4] =~ m/^\([ATCG]{3}\)$/i;

  # strip spaces from either end (the coordinates seem to have post-spaces)
  #map { s/^\s+// } @x;
  #map { s/\s+$// } @x;

  $count++;
  my $ID = make_id('trna', $count);

#  print STDERR "[$x[0]] Found $ID $x[1]$x[4] in $seqid\n" if $verbose;
  inform "[$count] Found $ID $x[1]$x[4] in $seqid" if $verbose;
  
  $x[2] =~ m/(c)?\[(\d+),(\d+)\]/;
  my $strand = defined $1 ? -1 : +1;
  my $start = $2;
  my $end = $3;

  my $f = Bio::SeqFeature::Generic->new( 
    -seq_id     => $seqid,
    -source_tag => $source,
    -primary    => 'tRNA',
    -start      => $start,
    -end        => $end,
    -strand     => $strand,
    -score      => undef,
    -frame      => 0,
    -tag        => {
      'ID' => $ID,
      'product' => "$x[1]$x[4]",
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
    {OPT=>"source=s", VAR=>\$source, DEFAULT=>'aragorn', DESC=>"tRNA prediction software"},
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
  print "Usage: $0 [options] < aragorn.out.file\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
