#!/usr/bin/env perl
use strict;

my(@Options, $verbose, $singles, $pairs);
setOptions();

$singles ||= '/dev/null';
open my $sing_io, '>', $singles;
$pairs ||= '/dev/null';
open my $pair_io, '>', ($pairs || '/dev/null');

my %seen;
my $nread = 0;
my $npass = 0;
my $loners = 0;

print STDERR "Reading: @ARGV\n";
print STDERR "Writing pairs: $pairs\n";
print STDERR "Writing singles: $singles\n";

# Correct orientaion is   =R3=>.............=F3=>
#                    or   <=F3=.............<=R3=
# so  R3 is /1 
# and F3 is /2 in velvet/nesoni speak

while (my $id = <ARGV>) {
#  print STDERR "$id";
  next unless $id =~ m/^>/;  
  $nread++;
  print STDERR "\rProcessed: $nread" if ($nread % 9871) == 0;
  my $seq = <ARGV> or last;
  if ($id =~ m/^>(\S+?)_([FR])\d$/) {
    my $stem = $1;  
    if (my $right = $seen{'F'}{$stem}) {
      print $pair_io ">$stem/1\n$seq>$stem/2\n$right";
      delete $seen{'F'}{$stem};
      $npass+=2;
    }
    elsif (my $left = $seen{'R'}{$stem}) {
      print $pair_io ">$stem/1\n$left>$stem/2\n$seq";
      delete $seen{'R'}{$stem};
      $npass+=2;
    }
    else {
      $seen{$2}{$1} = $seq;
    }
  }
  else {
    print $sing_io "$id\n$seq";
  }
}


for my $dir ('R', 'F') {
  for my $id (keys %{$seen{$dir}}) {
    print $sing_io ">${id}_$dir\n",$seen{$dir}{$id};
    $loners++;
  }
}

print STDERR "\rRead $nread, in pairs $npass, singles $loners\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbosity"},
    {OPT=>"s|singles=s", VAR=>\$singles, DEFAULT=>'/dev/null', DESC=>"Save singletons to this file"},
    {OPT=>"p|pairs=s",   VAR=>\$pairs, DEFAULT=>'/dev/null', DESC=>"Save pairs to this file"},
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
  print "Usage: $0 -s single.csfa -p pairs.csfa Reads_R3.csfa Reads_F3.csfa ...\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
