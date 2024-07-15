#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $chars, $flank);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
#my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $run=0;

#COUNT	SEQUENCE	POSITION	LEFT	LENGTH	CHARS	RIGHT
#1	VBC_Sal1296_C1	123758	actcagcgccgttgac	1312	N	ggcgccgctgacggtc
#2	VBC_Sal1296_C1	376279	agatcccggctaaatt	1678	N	ggaggttaatgctgag
#3	VBC_Sal1296_C1	469147	ttttgatatcgcctcg	63	N	gcgcgtctaccaattc
#4	VBC_Sal1296_C1	532826	gtcgcgataatggata	1101	N	ttcactattaactttg
#5	VBC_Sal1296_C1	600312	ccggatggcgacgcat	1229	N	gctcacggactgtgta
#6	VBC_Sal1296_C1	732355	ccacattatgcgcgac	1229	N	tttacgtcgtcatccg
#7	VBC_Sal1296_C1	777539	atggcggcggtggtgg	1	N	atatggacaacaccga

print join("\t",qw(COUNT SEQUENCE POSITION LEFT LENGTH CHARS RIGHT)),"\n";

while (my $seq = $in->next_seq) 
{
  my $s = $seq->seq;

  # pos($s) is the position (base-0) of the first non-char AFTER run (match)
  # but remember DNA is (base-1) ... sigh
  
  while ($s =~ m/([$chars]+)/gi) {
      print join("\t", 
      ++$run,
      $seq->id, 
      pos($s)-length($1)+1,  
      substr($s, pos($s)-length($1)-$flank, $flank),
      length($1),
      $chars,
      substr($s, pos($s), $flank),
    ),"\n";
  }

} 

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"c|chars=s",  VAR=>\$chars, DEFAULT=>'N', DESC=>"Characters to consider in runs (case-insenstive)"},
    {OPT=>"f|flank=i",  VAR=>\$flank, DEFAULT=>16, DESC=>"Flanking length to show"},
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
  print "Usage: $0 [options] long_scaffolds.fasta > shorter.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
