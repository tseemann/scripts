#!/usr/bin/env perl
use strict;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fasta write_fasta assert);

my(@Options, $verbose, $left, $right, $auto);
setOptions();

my $THROB = 7931;

print STDERR "Reading: @ARGV\n";

my $basename = $ARGV[0] || 'output';
if ($auto) {
  print STDERR "Using basename for --auto: $basename\n";
  $left = "$basename.left";
  $right = "$basename.right";
}
  
$left ||= '/dev/null';
print STDERR "Writing left: $left\n";
open my $left_io, '>', $left;

$right ||= '/dev/null';
print STDERR "Writing right: $right\n";
open my $right_io, '>', $right;

my $nread=0;
my $seq;

while ( not eof () ) {
  $seq = read_fasta(\*ARGV);
  write_fasta( (++$nread % 2 ? $left_io : $right_io) , $seq );
  print STDERR "\rProcessing: $nread" if $nread % $THROB == 0;
}

print STDERR "\rProcessed: $nread     \n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbosity"},
    {OPT=>"l|left=s", VAR=>\$left, DEFAULT=>'/dev/null', DESC=>"Left sequence file"},
    {OPT=>"r|right=s",   VAR=>\$right, DEFAULT=>'/dev/null', DESC=>"Right sequence file"},
    {OPT=>"a|auto!",   VAR=>\$auto, DEFAULT=>0, DESC=>"Auto: --left INPUT.left --right INPUT.right"},
  );

  #(@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [-l left.fa] [-r right.fa] interleaved.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
