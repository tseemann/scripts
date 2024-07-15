#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $uniqueids, $uppercase, $lowercase, $rename);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my %seen;

my $nread=0;
my $nwrote=0;
my $renamed=0;

while (my $seq = $in->next_seq) {
  my $id = $seq->display_id;
  $nread++;
  if ($seq->length <= 0) {
    print STDERR "Warning: >$id has no sequence, skipping.\n" if $verbose;
  }
  else {
    if ($uniqueids and $seen{$id}++) {
      $seq->display_id("$id.$renamed");
      $renamed++;
    }
    elsif ($rename) {
      $seq->display_id( sprintf "$rename", $nread );
      $renamed++;
    }
    $seq->seq( uc $seq->seq ) if $uppercase;
    $seq->seq( lc $seq->seq ) if $lowercase;
    $out->write_seq($seq);
    $nwrote++;
  }
} 
print STDERR "Read $nread sequences, wrote $nwrote, renamed $renamed, skipped ", $nread-$nwrote,".\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"uc!",  VAR=>\$uppercase, DEFAULT=>0, DESC=>"Uppercase sequences"},
    {OPT=>"lc!",  VAR=>\$lowercase, DEFAULT=>0, DESC=>"Lowercase sequences"},
    {OPT=>"u|uniqueids!",  VAR=>\$uniqueids, DEFAULT=>0, DESC=>"Modify IDs to ensure unique"},
    {OPT=>"r|rename=s",  VAR=>\$rename, DEFAULT=>'', DESC=>"Rename sequences using this template eg. 'contig%04d'"},
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
  print "Usage: $0 [options] dodgy.fasta > fixed.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
