#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;

my(@Options, $debug, $format, $gffver, $incseq, $idtag, $nametag);
setOptions();

my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);
print "##gff-version $gffver\n";

print STDERR "Loading: @ARGV ... (please be patient)\n";
my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>$format);
my @seq;
while (my $seq = $in->next_seq) {
  push @seq, $seq;
  printf "##sequence-region %s 1 %d\n", $seq->id, $seq->length;
}
printf STDERR "Loaded %d sequences\n", 0+@seq;

my $counter=0;
for my $seq (@seq) {
  for my $f ($seq->get_SeqFeatures) {
    $counter++;
    print STDERR "\rProcessing: ", $seq->display_id, " | $counter";
    my $id = $f->has_tag($idtag) ? ($f->get_tag_values($idtag))[0] : "SEQ$counter";
    $f->remove_tag('ID') if $f->has_tag('ID');
    $f->add_tag_value('ID', $id);

    my $name = $f->has_tag($nametag) ? ($f->get_tag_values($nametag))[0] : $id;
    $f->remove_tag('Name') if $f->has_tag('Name');
    $f->add_tag_value('Name', $name);
    
    print $f->gff_string($gff_factory),"\n";
  }
}
print STDERR "\rConverted $counter features into GFF $gffver\n";

if ($incseq) {
  print STDERR "Appending FASTA sequences\n";
  print "##FASTA\n";
  my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'fasta');
  for my $seq (@seq) {
    $out->write_seq($seq);
  }
}
else {
  print STDERR "Not adding sequences, specify --incseq for that.\n";
}
print STDERR "Done.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"format=s",  VAR=>\$format, DEFAULT=>'genbank', DESC=>"Input format"},
    {OPT=>"gffver=f",  VAR=>\$gffver, DEFAULT=>3, DESC=>"Output GFF format: 1 | 2 | 2.5 | 3"},
    {OPT=>"incseq!",  VAR=>\$incseq, DEFAULT=>0, DESC=>"Include sequence (GFF3 only)"},
    {OPT=>"idtag=s",  VAR=>\$idtag, DEFAULT=>'locus_tag', DESC=>"What to use as the GFF ID"},
    {OPT=>"nametag=s",  VAR=>\$nametag, DEFAULT=>'gene', DESC=>"What to use as the GFF Name"},
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
  print "Usage: $0 [options] < file.gbk > file.gff\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
