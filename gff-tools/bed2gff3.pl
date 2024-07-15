#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

my(@Options, $verbose, $source, $ftype, $fasta);
setOptions();

print "##gff-version 3\n";

my $seqdata;
open my $seqdata_fh, '>', \$seqdata;

if ($fasta) {
  my $in = Bio::SeqIO->new(-file=>$fasta, -format=>'fasta');
  my $out = Bio::SeqIO->new(-fh=>$seqdata_fh, -format=>'fasta');
  while (my $seq = $in->next_seq) {
    printf "##sequence-region %s 1 %d\n", $seq->id, $seq->length;
    $out->write_seq($seq);
  }
  print "###\n";
}

while (<>) {
  next if m/^#/;
  chomp;
  my @x = split /\t/;
  next if @x < 3;
  next unless $x[1] =~ m/^\d+$/;
  next unless $x[2] =~ m/^\d+$/;
  my($seqid, $start, $stop) = @x;
  $start = $start+1;
  my $strand = '+';
  if ($start > $stop) {
    ($start,$stop) = ($stop,$start);
    $strand = '-';
  }
  
  print join("\t",
    $seqid,
    $source,
    $ftype,
    $start,
    $stop,
    '.',
    '.',
    $strand,
    '.',
    ''
  ),"\n";
}  

if ($fasta) {
  print "##FASTA\n";
  print $seqdata;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"source=s", VAR=>\$source,  DEFAULT=>'bed2gff', DESC=>"What to put in column 2 of GFF"},
    {OPT=>"ftype=s", VAR=>\$ftype,  DEFAULT=>'misc_feature', DESC=>"What to put in column 3 of GFF"},
    {OPT=>"fasta=s", VAR=>\$fasta,  DEFAULT=>'', DESC=>"Include this FASTA file"},
#    {OPT=>"ids=s", VAR=>\$ids,  DEFAULT=>'Feature%d', DESC=>"Add ID tag with this format"},
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] file.bed > file.gff\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

