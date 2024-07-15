#!/usr/bin/env perl
use warnings;
use strict;
#use Regexp::Common;

#----------------------------------------------------------------------
my(@Options, $verbose);
setOptions();

# <operator>
my $op = shift @ARGV;
$op =~ m/^([\+\-\/\*])$/ or die "Invalid <operation> parameter";

# <listfiles>
open my $lfh, '<', (shift @ARGV);
open my $rfh, '<', (shift @ARGV);

while (1) {
  my $left = <$lfh>;
  last if not defined $left; 
  chomp $left;
  print STDERR "L=$left\n" if $verbose;
  my $right = <$rfh>;
  last if not defined $right;
  chomp $right;
  print STDERR "R=$right\n" if $verbose;
  my $result;
  eval "\$result = \$left $op \$right";
  print "$result\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
  );

  (@ARGV < 3) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage:\n  $0 [options] <operation> <listfile1> <listfile2>\n";
  print "  <operator> can be: + - * /\n";
  print "  <listfile> is a text file, one number per line.\n";
  print "Options:\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
