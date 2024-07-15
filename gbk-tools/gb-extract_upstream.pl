#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use List::Util qw(sum);

my(@Options, $verbose, $ftype, $nomask, $minlen, $gap, $len, $overlap, $inc_desc);
setOptions();

#tRNA            655..730
#                /gene="tRNA-Leu(UUR)"
#                /anticodon=(pos:678..680,aa:Leu)
#                /product="transfer RNA-Leu(UUR)"

my $gbk = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'genbank');
my $fna = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'fasta');
my $count=0;
my $num=0;

while (my $seq = $gbk->next_seq) {
  if (not $overlap) {
    mask_features($seq, $ftype, 'N')
  }
  my $L = $seq->length;
  for my $f (grep { $_->primary_tag eq $ftype } $seq->get_SeqFeatures) {
    my $id = TAG($f, 'locus_tag');
    print STDERR "\rProcessing: $id";
    $num++;
    my $desc = TAG($f, 'product') || TAG($f, 'gene');
    my $us;
    if ($f->strand < 0) {
      my($begin,$end) = ($f->end+1+$gap, $f->end+$len+$gap); 
      next if $begin < 1 or $end > $L;
      $us = $seq->trunc($begin, $end)->revcom;
    } 
    else {
      my($begin,$end) = ($f->start-$len-$gap, $f->start-1-$gap);
      next if $begin < 1 or $end > $L;
      $us = $seq->trunc($begin, $end);
    }
    my $dna = $us->seq;
    $dna =~ s/^.*N+//g;
    next unless $dna;              # skip if empty
    next if length $dna < $minlen; # skip if too short
    $us->seq($dna);
    $us->display_id(TAG($f,'locus_tag'));
    $desc .= " <upstream:".-($gap+1+$us->length)."..".-($gap+1).">";
    $desc .= " [".$seq->desc."]" if $inc_desc;
    $us->desc($desc);
    $fna->write_seq($us);
    $count++;
  }
}
print STDERR "\rWrote $count of $num upstream $ftype regions ~$len bp.\n";

#----------------------------------------------------------------------

sub mask_features {
  my($seq, $type, $char) = @_;
  my $dna = $seq->seq;
  for my $f (grep { $_->primary_tag eq $type } $seq->get_SeqFeatures) {
    substr $dna, $f->start, $f->length, ${char}x$f->length;
  }
  $seq->seq($dna);
}

#----------------------------------------------------------------------

sub TAG {
  my($f, $tag) = @_;
  return '' unless $f->has_tag($tag);
  return join(';', $f->get_tag_values($tag));
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"f|ftype=s",  VAR=>\$ftype, DEFAULT=>'CDS', DESC=>"Upstream of these"},
    {OPT=>"d|desc!",  VAR=>\$inc_desc, DEFAULT=>0, DESC=>"Include sequence description in FASTA desc"},
    {OPT=>"g|gap=i",  VAR=>\$gap, DEFAULT=>12, DESC=>"Don't include these immediate upstream bases"},
    {OPT=>"l|minlen=i",  VAR=>\$minlen, DEFAULT=>1, DESC=>"Minimum upstream sequence length"},
    {OPT=>"o|overlap!",  VAR=>\$overlap, DEFAULT=>0, DESC=>"Allow overlap into preceding feature"},
    {OPT=>"l|len=i",  VAR=>\$len, DEFAULT=>400, DESC=>"Maximum upstream sequence length"},
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
  print "Synopsis: Extracts upstream sequences of a given feature type\n";
  print "Usage: $0 [options] genome.gbk > upstream.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
