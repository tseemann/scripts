#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Data::Dumper;

use constant GFF_VERSION => 3;

my(@Options, $verbose, $pattern, $maxdist, $minlen, $feature, $idtag, $label, $incseqs, $incmainseq);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');

# pos($s) is the position (base-0) of the first non-char AFTER run (match)
# but remember DNA is (base-1) ... sigh

print "##gff-version ", GFF_VERSION, "\n";

my @seqnames;
my %inputseqs;

while (my $seq = $in->next_seq) 
{
  #load the sequences into a hash for later use if necessary
  if($incmainseq){
    push @seqnames, $seq->id;
    $inputseqs{$seq->id} = $seq;
  }

  # the DNA sequence string
  my $s = $seq->seq;
  
  # find all regions matching pattern and keep track in list
  my @r;
  while ($s =~ m/(${pattern}+)/g) {
#    push @r, [ pos($s), pos($s)+length($1)-1 ];
    push @r, [ pos($s)-length($1)+1, pos($s) ];
  }
  printf STDERR "%s: found %d regions (pattern=$pattern)\n", 
    $seq->display_id, scalar(@r), $pattern;    
  print STDERR Dumper(\@r) if $verbose;
  next unless @r;
  
  # merge adjoining list items that are within $maxdist
  my @p;
  my $p = shift @r;
  for my $r (@r) {
    if ($r->[0] - $p->[1] <= $maxdist) {
      #$p = [ $p->[0], $r->[1] ];
      $p->[1] = $r->[1]; # extend end across
    }
    else {
      push @p, $p;
      $p = [ @$r ];
    }
  }
  push @p, $p if $p; # leftover one
  printf STDERR "%s: merged into %d regions (mergedist=%d)\n", 
    $seq->display_id, scalar(@p), $maxdist;    
  print STDERR Dumper(\@p) if $verbose;
  next unless @p;
  
  # eliminate regions smaller than $minlen
  
  my @q = grep { $_->[1] - $_->[0] >= $minlen } @p;
  printf STDERR "%s: reduced to %d regions (minlen=%d)\n", 
    $seq->display_id, scalar(@q), $minlen;    
  print STDERR Dumper(\@q) if $verbose;
  next unless @q;

#  print Dumper(\@q);

  # output as GFF
  my $count=0;
  my $gff_factory = Bio::Tools::GFF->new(-gff_version=>GFF_VERSION);
  for my $q (sort { $a->[0] <=> $b->[0] } @q) {
    my $id = sprintf "$idtag%05d", ++$count;
    my $feat = Bio::SeqFeature::Generic->new(
      -seq_id  => $seq->display_id,
      -source  => 'vbc',
      -primary => $feature,
      -start   => $q->[0],
      -end     => $q->[1],
      -strand  => '.',
#      -score   => '.',
      -tag => {
        ID      => $id,
	product => "$label - ".abs($q->[1]-$q->[0]+1)." bp",
	colour  => '200 0 0',
      }
    );
    print $feat->gff_string($gff_factory),"\n";
    $seq->add_SeqFeature($feat);
  }
  
  if ($incseqs) {
    my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'fasta');
    print "##FASTA\n";
    for my $f ($seq->get_SeqFeatures) {
      my $s = $f->seq;
      $s->display_id( $f->get_tag_values('ID') );
      $s->desc( $f->get_tag_values('product') );
      $out->write_seq($s);
    }
  }

  if($incmainseq){
    my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'fasta');
    unless($incseqs){ print "##FASTA\n"; }
    foreach my $id (@seqnames){
      $out->write_seq($inputseqs{$id});
    }
  }
    
  printf STDERR "%s: wrote %d %s %s features\n", 
    $seq->display_id, scalar(@q), $feature, $label;
  
} 

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"p|pattern=s", VAR=>\$pattern,DEFAULT=>'[atgcn]',  DESC=>"Pattern to count as regions (case sensitive)"},
    {OPT=>"d|mergedist=i", VAR=>\$maxdist,DEFAULT=>25,  DESC=>"Join regions this close together"},
    {OPT=>"m|minlen=i",    VAR=>\$minlen, DEFAULT=>250, DESC=>"Minimum region size to allow"},
    {OPT=>"f|feature=s", VAR=>\$feature, DEFAULT=>'misc_feature',  DESC=>"Feature type for GFF output"},
    {OPT=>"t|idtag=s", VAR=>\$idtag, DEFAULT=>'noncore',  DESC=>"Prefix for feature IDs"},
    {OPT=>"l|label=s", VAR=>\$label, DEFAULT=>'non-core genome',  DESC=>"Label for /product of each feature"},
    {OPT=>"s|sequences!", VAR=>\$incseqs, DEFAULT=>0,  DESC=>"Include sequences in FASTA block at end of GFF3 output"},
    {OPT=>"i|include!", VAR=>\$incmainseq, DEFAULT=>0,  DESC=>"Include the original sequence in FASTA block at end of GFF3 output"},
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
  print "Usage: $0 [options] file.fasta > file.pattern.gff\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

