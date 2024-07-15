#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;
use List::Util qw(sum);

my(@Options, $debug, $mincov);
setOptions();

my %len_of;
my %depth_of;
my %low;

for my $bamfn (@ARGV) {
  -r $bamfn or die "can't read bam file: $bamfn";
  foreach (qx(samtools view -H \Q$bamfn\E)) {
    next unless m/^\@SQ\s+SN:(\S+)\s+LN:(\d+)$/;
    $len_of{$1} = $2;
    $low{$1} = 0;
    $depth_of{$1} = 0;
  }
  my $totbp = sum(values %len_of);
  open DEPTH, "samtools depth \Q$bamfn\E |";
  while (<DEPTH>) {
    my @x = split m/\t/;
    $depth_of{$x[0]} += $x[2];
    $low{$x[0]}++ if $x[2] <= $mincov;
    print STDERR "\rProcessing: $./$totbp" if $. % 9871 == 0;
  }              
  print STDERR "\n";
  print "Sequence\tLength\tAvgDepth\tDepths<=$mincov\n";
  for my $id (keys %depth_of) {
    printf "%s\t%d\t%.1f\t%d\n", 
     $id, 
     $len_of{$id}, 
     $depth_of{$id}/$len_of{$id}, 
     $low{$id},
  }
}

#print Dumper(\%len_of, \%depth_of);


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"mincov=i",  VAR=>\$mincov, DEFAULT=>0, DESC=>"Min coverage to count as uncovered"},
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
  print "Usage: $0 [options] <file.bam>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
