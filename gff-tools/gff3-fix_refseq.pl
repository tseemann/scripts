#!/usr/bin/env perl
use strict;
use warnings;

my(@Options, $verbose, $gffver);
setOptions();

my @line = <ARGV>;
my @head = grep {  m/^##[^#]/ } @line;
my @feat = grep { !m/^#/ } @line;

printf STDERR "Read %d lines, %d header, %d features.\n", 0+@line, 0+@head, 0+@feat;

#gff-version 3
##sequence-region NC_007790.1 1 3125
##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=451515
##sequence-region NC_007791.1 1 4439
##sequence-region NC_007792.1 1 37136

my %seen;
foreach (@head) {
  print if not $seen{$_}++;
}

#NC_007790.1 RefSeq gene 2093 2674 . - . ID=gene3;Name=SAUSA300_pUSA010004;Dbxref=GeneID:3912720;gbkey=Gene;locus_tag=SAUSA300_pUSA010004
#NC_007790.1 RefSeq CDS 2093 2674 . - 0 ID=cds3;Name=YP_492682.1;Parent=gene3;Dbxref=Genbank:YP_492682.1,GeneID:3912720;gbkey=CDS;product=hypothetical protein;protein_id=YP_492682.1;transl_table=11

print "###\n";
foreach (@feat) {
  my @x = split m/\t/;
  next if $x[2] =~ m/^(CDS|exon)$/;
  $x[8] =~ s/ID=\w+;//;
  $x[8] =~ s/locus_tag=/ID=/;
  $x[8] =~ s/gbkey=\w+;//;
  print join("\t", @x);
}

print "###\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
  print "Usage: $0 [options] bacteria_refseq.gff ... > proper_ids.gff\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

