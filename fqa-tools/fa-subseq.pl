#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use File::Basename;

my(@Options, $debug, $begin, $end, $rename);
setOptions();

my $in  = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

while (my $seq = $in->next_seq) {
    my $p1 = $begin ? $begin : 1;
    my $p2 = $end   ? $end   : $seq->length;
    print STDERR $seq->display_name, " $p1..$p2\n";
    $seq =  $seq->trunc($p1, $p2);
    $seq->id($seq->id.":$p1-$p2") if $rename;
    $out->write_seq($seq);
} 

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"begin|start=i",  VAR=>\$begin, DEFAULT=>0, DESC=>"Start (0=existing start)"},
    {OPT=>"end=i",  VAR=>\$end, DEFAULT=>0, DESC=>"End (0=existing end)"},
    {OPT=>"rename!",  VAR=>\$rename, DEFAULT=>0, DESC=>"Append :being-end to ID"},
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
  my $exe = basename($0);
  print "Usage: $exe [options] in.fasta > out.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
