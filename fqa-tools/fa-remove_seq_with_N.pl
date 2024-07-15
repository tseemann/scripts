#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $ids);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $nread=0;
my $nwrote=0;

while (my $seq = $in->next_seq) {
  $nread++;
#  print STDERR "Replacing Ns in: ",$seq->display_id,"\n" if $verbose;
  next if $seq->seq =~ m/N/i;
  if ($ids) {
    print $seq->id,"\n";
  }
  else {
    $out->write_seq($seq);
  }
  $nwrote++;
} 
printf STDERR "Read $nread sequences, wrote $nwrote (discarded %d)\n", ($nread-$nwrote);

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"ids!",  VAR=>\$ids, DEFAULT=>0, DESC=>"Output IDs rather than FASTA"},
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
