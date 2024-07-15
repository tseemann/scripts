#!/usr/bin/env perl
use strict;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fasta write_fasta assert);

my(@Options, $verbose, $singles, $pairs, $auto, $newbler);
setOptions();

my $basename = $ARGV[0] || 'output';
if ($auto) {
  print STDERR "Using basename for --auto: $basename\n";
  $singles = "$basename.singles";
  $pairs = "$basename.pairs";
}

$singles ||= '/dev/null';
open my $sing_io, '>', $singles;
$pairs ||= '/dev/null';
open my $pair_io, '>', $pairs;
#open my $pair_io, '>', ($pairs || '/dev/null');

my %seen;
my $nread = 0;
my $npass = 0;
my $loners = 0;
my $lib = 'Pairs';

print STDERR "Reading: @ARGV\n";
print STDERR "Writing pairs: $pairs\n";
print STDERR "Writing singles: $singles\n";

while ( not eof () ) {   # perldoc -f eof
  my $seq = read_fasta(\*ARGV);
  $nread++;
  print STDERR "\rProcessed: $nread" if $nread % 97131 == 0;
#  print "seq[0]=|$seq->[0]|\n" if $seq->[0] =~ m{/};
  if ($seq->[0] =~ m{^(\S+?)/([12])}) {
    my $stem = $1;  
    # template=ZV76N:99:946 dir=F library=PGM
    my $Ldesc = $newbler ? " template=$stem dir=F library=$lib" : "";
    my $Rdesc = $newbler ? " template=$stem dir=R library=$lib" : "";
    if (my $left = $seen{'1'}{$stem}) {
      print $pair_io ">$stem/1$Ldesc\n$left\n>$stem/2$Rdesc\n",$seq->[1],"\n";
      delete $seen{'1'}{$stem};
      $npass+=2;
    }
    elsif (my $right = $seen{'2'}{$stem}) {
      print $pair_io ">$stem/1$Rdesc\n",$seq->[1],"\n>$stem/2$Rdesc\n$right\n";
      delete $seen{'2'}{$stem};
      $npass+=2;
    }
    else {
      $seen{$2}{$1} = $seq->[1];
      $loners++;
    }
  }
  else {
    print $sing_io ">",$seq->[0],"\n",$seq->[1],"\n";
  }
}
print STDERR "\n";

for my $dir (1..2) {
  for my $id (keys %{$seen{$dir}}) {
    print $sing_io ">$id\n",$seen{$dir}{$id},"\n";
    $loners++;
  }
}

print STDERR "Read $nread, in pairs $npass, singles $loners\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbosity"},
    {OPT=>"singles=s",VAR=>\$singles, DEFAULT=>'/dev/null', DESC=>"Save singletons to this file"},
    {OPT=>"pairs=s",  VAR=>\$pairs, DEFAULT=>'/dev/null', DESC=>"Save pairs to this file"},
    {OPT=>"auto!",    VAR=>\$auto, DEFAULT=>0, DESC=>"Auto: --singles INPUT.singles --pairs INPUT.pairs"},
    {OPT=>"newbler!", VAR=>\$newbler, DEFAULT=>0, DESC=>"Add Newbler meta-data to fasta description for pairs"},
  );

  (@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 -s unmapped.single.fa -p unmapped.pair.fa unmapped.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
