#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $replace, $list);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $nread=0;
my $oldbp=0;
my $newbp=0;
my $runs=0;

while (my $seq = $in->next_seq) {
  print STDERR "Replacing Ns in: ",$seq->display_id,"\n" if $verbose;
  my $s = $seq->seq;
  $oldbp += length($s);
  my $n = ($s =~ s/N+/$replace/gi);
  $runs += $n;
#  print STDERR "n=$n\n";
  $newbp += length($s);
  if ($list) {
    print join("\t",$seq->id,$seq->length,$n),"\n";
  }
  else {
    $seq->seq($s);
    $out->write_seq($seq);
  }
  $nread++;
} 
print STDERR "Read $nread sequences, replaced $runs N-runs with '$replace'. Had $oldbp bp, now $newbp. Removed ",($oldbp-$newbp)," Ns\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"r|replace=s",  VAR=>\$replace, DEFAULT=>'N', DESC=>"Replace runs of N with this"},
    {OPT=>"l|list!",  VAR=>\$list, DEFAULT=>0, DESC=>"List N-runs rather than print FASTA"},
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
