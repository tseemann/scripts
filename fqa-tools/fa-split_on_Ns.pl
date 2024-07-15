#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $sep, $minlen);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $nread=0;
my $nwrote=0;

while (my $seq = $in->next_seq) {
  print STDERR $seq->display_id,"\n" if $verbose;
  $nread++;
  if ($seq->seq =~ m/N/i) {
    my $count=0;
    for my $s (split m/N+/i, $seq->seq) {
      next unless length($s) >= $minlen;
      $count++;
      print STDERR "\t$count\n" if $verbose;
      $out->write_seq( Bio::Seq->new(
        -id=>$seq->display_id.$sep.$count, 
	-seq=>$s,
#	-desc=>$seq->desc
      ));	
      $nwrote++;
    }
  }
  else {
    $out->write_seq($seq);
    $nwrote++;
  }  
}

print STDERR "Read $nread, wrote $nwrote sequences.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"sep=s", VAR=>\$sep, DEFAULT=>'_', DESC=>"Separator string for IDs"},
    {OPT=>"l|minlen=i", VAR=>\$minlen, DEFAULT=>31, DESC=>"Minimum length"},
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
  print "Usage: $0 [options] scaffolds.fasta > contigs.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
