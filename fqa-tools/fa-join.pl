#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;

my(@Options, $verbose, $sep, $gffout, $noends);
setOptions();

my $in  = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $gff_factory = Bio::Tools::GFF->new(-gff_version=>3);

open my $gff, '>', $gffout;

my $lsep = length($sep);

my($start,$cat) = $noends ? (1,'') : (1+$lsep,$sep);
my $ncat = 0;
my $previd = '';

while (my $seq = $in->next_seq) 
{
#  next unless $seq->alphabet eq 'dna';

  my $end = $start + $seq->length - 1;

  my $f = Bio::SeqFeature::Generic->new(
    -start       => $start,
    -end         => $end,
    -seq_id      => $seq->display_name,
    -source_tag  => $0,
    -primary_tag => 'source',
    -tag         => {
      ID      => $seq->display_name,
      product => $seq->desc,
#      note    => "seq_no_".($ncat+1),
    }
  );
  print $gff $f->gff_string($gff_factory),"\n";
 
  $cat .= $seq->seq;
  $cat .= $sep;
  $start = $end + $lsep + 1;
  $ncat++;
  $previd = $seq->display_name;
} 
close $gff;

$cat .= $sep unless $noends;

my $new = Bio::PrimarySeq->new(
  -id   => "join_${ncat}_seq",
  -desc => "Joined $ncat sequences with $sep separator, ".length($cat)." bp.",
  -seq  => $cat,
);

my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
$out->write_seq($new);

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"verbose info"},
    {OPT=>"sep=s",  VAR=>\$sep, DEFAULT=>'NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN', DESC=>"Separator sequence"},
    {OPT=>"gffout=s",  VAR=>\$gffout, DEFAULT=>'/dev/null', DESC=>"Write GFF for contig coordinates joins"},
    {OPT=>"noends!",  VAR=>\$noends, DEFAULT=>0, DESC=>"Don't put joiner at start and end"},
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
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
