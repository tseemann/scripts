#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqFeature::Generic;
use List::Util qw(min max);

my(@Options, $verbose, $gff_version, $source);
setOptions();

my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gff_version);

# orf00002     1390      356  -2     9.64

my $seqid='';
my $count=0;

while (<ARGV>) 
{
  chomp;
  my @x = split ' ';
  
  if (m/^>(\S+)/) {
    $seqid = $1;
    next;
  }
  
  next unless @x == 5 and $x[1] =~ m/^\d+$/ and $x[2] =~ m/^\d+$/;

  my $ID = sprintf("cds%05d", ++$count); 
  print STDERR "[$count] Found $ID in $x[0]\n" if $verbose;

  my($strand,$left,$right) = $x[1] <= $x[2] ? (+1,$x[1],$x[2]) : (-1,$x[2],$x[1]);
  my $partial=0;
  if ($left < 1) { $left+=3; $partial=1 }

#  can't know the end of the sequence...
#  if ($right > $seq[$N]->length) { $right-=3; $partial=1 }

  my $f = Bio::SeqFeature::Generic->new( 
    -seq_id     => $seqid,
    -source_tag => $source,
    -primary    => 'CDS',
    -start      => $left,
    -end        => $right,
    -strand     => $strand,
#    -score      => $x[4],
    -tag        => {
      'ID' => $ID,
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
    {OPT=>"source=s", VAR=>\$source, DEFAULT=>'Glimmer3', DESC=>"Glimmer prediction software"},
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
