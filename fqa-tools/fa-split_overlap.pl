#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use List::Util qw(min max);

my(@Options, $verbose, $length, $overlap, $phrap, $insert);
setOptions();

if (not defined $overlap or $overlap < 0) {
  $overlap = int($length/2);
}

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $delta = $length - $overlap;

my $count=0;
my $contig=0;
while (my $seq = $in->next_seq) {
  $contig++;
  print STDERR $seq->display_id,"\n" if $verbose;

  for ( my $start = 1 ; $start <= $seq->length ; $start += $delta ) {
    my $end = $start+$length-1;
    while ($end < $seq->length and $seq->subseq($start,$end) =~ m/N$/) {
      $end += $overlap;
      $start += $overlap;
    }
    $end = $seq->length if $end+$delta > $seq->length;   

    $end = min($seq->length, $end);
    
    my $read = $seq->trunc($start, $end);
    $count++;
    $read->display_id("$count.b");
    $read->desc($seq->display_id.":$start-$end");
    $out->write_seq($read);
  }
}

#print STDERR "Read $nread, wrote $nwrote sequences.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"l|length=i", VAR=>\$length, DEFAULT=>'3000', DESC=>"Sequence length"},
    {OPT=>"o|overlap=i", VAR=>\$overlap, DEFAULT=>'', DESC=>"Overlap length. -1 = half the length"},
#    {OPT=>"phrap!", VAR=>\$phrap, DEFAULT=>0, DESC=>"Use phrap IDs"},
#    {OPT=>"i|insert", VAR=>\$insert, DEFAULT=>0, DESC=>"Insert size (0=no pairing)"},    
  );

#  (@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 long.fasta > short.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
