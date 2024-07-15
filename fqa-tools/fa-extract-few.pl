#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $inverse);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*STDIN, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
my $nread=0;
my $nwrote=0;

my $pattern = join('|', @ARGV);

while (my $seq = $in->next_seq) {
  $nread++;
  my $match = ($seq->description =~ m/($pattern)/ or $seq->display_id =~ m/($pattern)/);
  print STDERR "Found match: ",$seq->display_id, " ", $seq->description, "\n" if $verbose;
  if ($match ^ $inverse) {  # rare use for XOR !
    $out->write_seq($seq);
    $nwrote++;
  }
} 

print STDERR "Read $nread sequences, wrote $nwrote, with pattern: $pattern\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"h|help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"v|inverse!",  VAR=>\$inverse, DEFAULT=>0, DESC=>"Output NON-matching sequences instead"},
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
  print "Usage: $0 [options] id1 [id2 ...] < input.fasta > output.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
