#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $unique);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
#my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my %seq_of;
my $nread=0;

while (my $seq = $in->next_seq) {
  print STDERR $seq->id, "\t", $seq->length, "\t", $seq->alphabet, "\n" if $verbose;
  $seq_of{$seq->id} .= $seq->seq;
  $nread++;
} 

my $nwrote=0;

for my $id (sort keys %seq_of) {
  print ">$id\n",$seq_of{$id},"\n";
  $nwrote++;
}

print STDERR "Read $nread, wrote $nwrote.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose info"},
#    {OPT=>"u|unique!",  VAR=>\$unique, DEFAULT=>0, DESC=>"Don't even output one representative for duplicates"},
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
  print "Usage: $0 [options] seq.fasta [ more.fasta ... ] > concatonid.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
