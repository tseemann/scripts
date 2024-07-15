#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Bio::Tools::IUPAC;
use Bio::Tools::SeqPattern;

my(@Options, $verbose, $inverse, $fasta, $countem, $revcom);
setOptions();

my $pattern = shift @ARGV;
print STDERR "Pattern: $pattern\n" if $verbose;
my $patt = Bio::Tools::SeqPattern->new(-SEQ=>$pattern, -TYPE=>'dna');
die "bad pattern '$pattern'" unless $patt->alphabet_ok;

my $regexp = $patt->expand;
$regexp .= '|'.$patt->revcom->expand if $revcom;
print STDERR "Regexp: $regexp\n" if $verbose;

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $num=0;
while (my $seq = $in->next_seq) {
  my $dna = $seq->seq;
  while ($dna =~ m/($regexp)/gi) {
    my $match = $1;
    my $pos = pos($dna);
    $num++;
    next if $countem;
    if ($fasta) {
      print ">", $seq->display_id," ",$seq->desc,"\n$match\n";
    }
    else {
      if ($seq->desc =~ /<upstream:(-\d+)\.\./) {
        # special coordinate transform - HACK!
        $pos = $pos + $1;
      }
      print join("\t",
        $seq->display_id, $match, $pos, $seq->desc
      ),"\n";
    }
  } 
} 
print "$num\n" if $countem;

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,  DESC=>"This help"},
    {OPT=>"d|debug!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
#    {OPT=>"v|inverse!",  VAR=>\$inverse, DEFAULT=>0, DESC=>"Return non-matches instead of matches"},
    {OPT=>"f|fasta!",  VAR=>\$fasta, DEFAULT=>0, DESC=>"Output in fasta"},
    {OPT=>"c|count!",  VAR=>\$countem, DEFAULT=>0, DESC=>"Just count"},
    {OPT=>"r|rc|revcom!",  VAR=>\$revcom, DEFAULT=>0, DESC=>"Check revcom of pattern too"},
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
  print "Usage: $0 [options] 'pattern' input.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
