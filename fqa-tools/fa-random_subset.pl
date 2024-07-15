#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $num, $prop);
setOptions();

die "please specify -n or -p" unless ($num or $prop);
die "bad -p $prop : must be between 0 and 1" if $prop and $prop < 0 or $prop > 1;

print STDERR "Counting no. of sequences in input...";
my($nseq) = qx(grep -c '>' @ARGV);
chomp $nseq;
$num = int($prop*$nseq) if $prop;
print STDERR " found $nseq, choosing random $num.\n";

my @pick = map { 0 } (1 .. $nseq);
my $picked=0;
while ($picked < $num) {
  my $i = int rand($nseq-1);
  next if $pick[$i];
  $pick[$i] = 1;
#  print STDERR "picked $i\n";
  $picked++;
}

my $in  = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
my $found=0;
#print STDERR "nseq $nseq\n";
for my $i (0 .. $nseq-1) {
#  print STDERR "seq $i found $found\n";
  my $seq = $in->next_seq;
  next unless $pick[$i];
  $out->write_seq($seq);
  $found++;
  last if $found >= $num;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"verbose info"},
    {OPT=>"n|num=i",  VAR=>\$num, DEFAULT=>0, DESC=>'Choose N sequences from input'},
    {OPT=>"p|proportion=f",  VAR=>\$prop, DEFAULT=>0, DESC=>'Choose this proportion of input sequences'},
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
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
