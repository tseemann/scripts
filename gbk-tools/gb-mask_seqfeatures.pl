#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

my(@Options, $verbose, $ftype, $informat, $outformat, $mask);
setOptions();

my $in  = Bio::SeqIO->new(-fh=>\*ARGV,   -format=>$informat);
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>$outformat);

$mask ||= 'N';
$mask = substr $mask, 0, 1;
print STDERR "Using mask: $mask\n";

my $count = 0;

my $inputs = @ARGV ? "@ARGV" : "(STDIN)";
print STDERR "Reading: $inputs\n";
while (my $seq = $in->next_seq) {
  $count++;
  print STDERR "[$count] ", $seq->id, " ... ";
  my $n = mask_features($seq, $ftype, $mask);
  $out->write_seq($seq);
  print STDERR "$n $ftype features masked.\n";
}
#print STDERR "\rWrote $count of $num upstream $ftype regions ~$len bp.\n";

#----------------------------------------------------------------------

sub mask_features {
  my($seq, $type, $char) = @_;
  my $dna = $seq->seq;
  my $count = 0;
  for my $f (grep { $_->primary_tag eq $type } $seq->get_SeqFeatures) {
    substr $dna, $f->start, $f->length, ${char}x$f->length;
    $count++;
  }
  $seq->seq($dna);
  return $count;
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
    {OPT=>"if|informat=s",  VAR=>\$informat, DEFAULT=>'genbank', DESC=>"Input format"},
    {OPT=>"of|outformat=s",  VAR=>\$outformat, DEFAULT=>'genbank', DESC=>"Output format"},
    {OPT=>"f|ftype=s",  VAR=>\$ftype, DEFAULT=>'rRNA', DESC=>"Mask these"},
    {OPT=>"m|maskchar=s",  VAR=>\$mask, DEFAULT=>'N', DESC=>"Mask character"},
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
  print "Synopsis: Masks sequence of a feature type with Ns\n";
  print "Usage: $0 [options] original.gbk > masked.gbk\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
