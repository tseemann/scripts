#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::FeatureIO;

my(@Options, $verbose, $ftype, $translate, $id_field, $desc_field, $suffix);
setOptions();

#----------------------------------------------------------------------

my %cds;
my %seq;
my $num=0;

my $in = Bio::FeatureIO->new(-format=>'gff', -version=>3, -fh=>\*ARGV);

while (my $f = $in->next_feature) {
  next unless $f->primary_tag eq $ftype;
  next unless $f->has_tag($id_field);
  my($ID) = $f->get_tag_values($id_field);
  push @{$cds{$ID}}, $f;
  print STDERR "\rProcessing: ", ++$num;
}

while (my $seq = $in->next_seq) {
  $seq{ $seq->display_id } = $seq;
  print STDERR "\rProcessing: ", ++$num;
}

print STDERR "\n";
my $out = Bio::SeqIO->new(-format=>'fasta', -fh=>\*STDOUT);

for my $ID (keys %cds) {
  my @cds = sort { $a->start <=> $b->start } @{$cds{$ID}};
  @cds = reverse @cds if $cds[0]->strand < 0;
  my $dna = '';
  for my $c (@cds) {
    $c->attach_seq( $seq{ $c->seq_id } );
    $dna .= $c->seq->seq;
  }
  my $pep = Bio::Seq->new(-id=>"$ID$suffix", -seq=>$dna);
  $pep = $pep->translate if $translate;
  $out->write_seq($pep);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"id=s", VAR=>\$id_field, DEFAULT=>'ID', DESC=>'GFF3 attribute to join on'},
    {OPT=>"suffix=s", VAR=>\$suffix, DEFAULT=>'', DESC=>'Append this to ID on output sequences'},
    {OPT=>"ftype=s", VAR=>\$ftype, DEFAULT=>'CDS', DESC=>"Feature to splice"},
    {OPT=>"translate!", VAR=>\$translate, DEFAULT=>0, DESC=>"Translate spliced feature?"},
#    {OPT=>"desc=s", VAR=>\$desc_field, DEFAULT=>'product', DESC=>'GFF3 field as Description'},
  );

#  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] file.gff3 > file.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

