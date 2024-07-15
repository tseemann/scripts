#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;

use FindBin;
use lib "$FindBin::Bin/";
use Prokka qw(inform read_gff3);

my(@Options, $verbose, $tags, $id, $desc, $prot, $trunc, $stop, $nometh);
setOptions();

my $gff = read_gff3(@ARGV);
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'fasta');

if ($tags =~ m/^(|source|all)$/i) {
  # just want all the sequences present in ##FASTA section
  $out->write_seq(values %{$gff->{'sequences'}});
}
else {
  for my $f (@{$gff->{'features'}}) 
  {
    if ($f->primary_tag =~ m/$tags/i) {
      inform $f->gff_string if $verbose;
      
      print Dumper($f); exit;
      
      my $seq = $f->spliced_seq(-phase => $f->phase);
#      $seq = $f->cds;
      
      if ($prot) {
	$seq = $seq->translate;
	$seq = $seq->trunc(1, $seq->length-1) if !$stop and $seq->length > 1 and $seq->seq =~ m/\*$/;
	if (!$nometh) {
          my $s = $seq->seq;
	  substr($s,0,1,'M');
	  $seq->seq($s);
	}
      }  

      if ($trunc and $trunc < $seq->length) {
	$seq = $trunc > 0 ? $seq->trunc(1, $trunc)
                          : $seq->trunc($seq->length+$trunc, $seq->length);
      }
      if ($id and $f->has_tag($id)) {
	$seq->display_id( ($f->get_tag_values($id))[0] );
      }
      if ($desc and $f->has_tag($desc)) {
	$seq->desc( ($f->get_tag_values($desc))[0] );
      }
      $out->write_seq($seq);
    }
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"t|tags=s", VAR=>\$tags,
        DEFAULT=>'^CDS$', DESC=>'Regexp of feature types to take part (Use "source" for whole genomic sequence)'},
    {OPT=>"id=s", VAR=>\$id, DEFAULT=>'ID', DESC=>'GFF3 field as >ID'},
    {OPT=>"desc=s", VAR=>\$desc, DEFAULT=>'product', DESC=>'GFF3 field as Description'},
    {OPT=>"p|prot!", VAR=>\$prot, DEFAULT=>0, DESC=>'Translate to protein'},
    {OPT=>"trunc=i", VAR=>\$trunc, DEFAULT=>0, DESC=>'Truncate to this size (-ve from end)'},
    {OPT=>"stop!", VAR=>\$stop, DEFAULT=>0, DESC=>"Leave STOP codon from end"},
    {OPT=>"nometh!", VAR=>\$nometh, DEFAULT=>0, DESC=>"Don't change START to M (Met)"},

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

