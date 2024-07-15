#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $id_re, $desc_re, $seq_re, $inverse, $noseq, $and_logic);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

while (my $seq = $in->next_seq) 
{
  my $id_ok = $id_re && $seq->display_id =~ m/$id_re/i;
  my $desc_ok = $desc_re && $seq->desc =~ m/$desc_re/i;
  my $seq_ok = $seq_re && $seq->seq =~ m/$seq_re/i;
  
  my $matches = (!$and_logic) ? ($id_ok || $desc_ok || $seq_ok) : ($id_ok && $desc_ok && $seq_ok);
#  my $matches = ($id_ok || $desc_ok || $seq_ok) : ($id_ok && $desc_ok && $seq_ok);
  $matches = !$matches if $inverse;
  
  if ($matches) {
    if ($noseq) {
      printf "%s %s\n", $seq->display_id, $seq->desc;
    }
    else {
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
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"i|id=s",  VAR=>\$id_re, DEFAULT=>'', DESC=>"ID regexp"},
    {OPT=>"d|desc=s",  VAR=>\$desc_re, DEFAULT=>'', DESC=>"Description regexp"},
    {OPT=>"s|seq=s",  VAR=>\$seq_re, DEFAULT=>'', DESC=>"Sequence regexp"},
    {OPT=>"v|inverse!",  VAR=>\$inverse, DEFAULT=>0, DESC=>"Return non-matches instead of matches"},
    {OPT=>"n|noseq!",  VAR=>\$noseq, DEFAULT=>0, DESC=>"Don't output full sequences, just ID & desc."},
    {OPT=>"a|and!",  VAR=>\$and_logic, DEFAULT=>0, DESC=>"Use AND rather than OR logic"},
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
  print "Usage: $0 [options]    input.fasta > matches.fasta\n";
  print "       $0 [options] -v input.fasta > non-matches.fasta\n";
  print "       $0 [options] -n input.fasta > matches.ids.txt\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
