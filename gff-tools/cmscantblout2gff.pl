#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
#use Bio::SearchIO;
use Data::Dumper;
use List::Util qw(min max);

##target name         accession query name                          accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
##------------------- --------- ----------------------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
#tRNA                 RF00005   NODE_148_length_2759_cov_152.375854 -          cm        1       71     1776     1704      -    no    1 0.32   0.0   61.8   1.9e-12 !   tRNA
#tRNA                 RF00005   NODE_148_length_2759_cov_152.375854 -          cm        1       71     1606     1532      -    no    1 0.33   0.0   55.7   8.9e-11 !   tRNA


my(@Options, $verbose, $gff_version, $source, $skip);
setOptions();

my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gff_version);
my $count=0;

while (<ARGV>) 
{
  if ($count==0 and !m/^\#target/) {
    die "Input does not look like Infernal 1.1 'cmscan --tblout' format";
  }
  $count++;
  my @x = split ' ';   # yes, split on _whitespace_, these are not tab separated!
#  print Dumper(\@x);
  next if $x[0] =~ m/^#/;
  next if $x[0] =~ m/$skip/;
  next unless defined $x[9];
  my $desc = join(' ', @x[17 .. $#x]) if @x > 16;
  my $f = Bio::SeqFeature::Generic->new( 
    -seq_id     => $x[2],
    -source_tag => 'Infernal',
    -primary    => 'ncRNA',
    -start      => min($x[7], $x[8]),
    -end        => max($x[7], $x[8]),
    -strand     => ($x[9] eq '+' ? +1 : -1),
    -score      => $x[15],
    -frame      => 0,
    -tag        => {
#        'ID' => $ID,
      'product' => "Rfam:$x[1] $desc",
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
    {OPT=>"skip=s", VAR=>\$skip, DEFAULT=>'tRNA|_rRNA_', DESC=>"Skip IDs with this pattern"},
    
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
