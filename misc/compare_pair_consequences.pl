#!/usr/bin/env perl
use strict;
use Data::Dumper;


my(@Options, $verbose);
setOptions();

my $fx = shift @ARGV;
my $fy = shift @ARGV;
my $x = load_nesoni_csq($fx);
my $y = load_nesoni_csq($fy);

#print Dumper($x);

my @kx = keys %$x;
my @ky = keys %$y;
delete @{$x}{@ky};
delete @{$y}{@kx};

print "*** $fx VERSUS $fy ***\n";
print "Changes in $fx (not in $fy)\n";
#print Dumper($x);
#print "1/2\n";
print values %$x;
#print "2/2\n";
print "Changes in $fy (not in $fx)\n";
#print Dumper($y);
print values %$y;

sub load_nesoni_csq {
  my $fname = shift;
  open my $fh, '<', $fname;
  my %h;
  while (<$fh>) {
#    chomp;
    my @f = split m/\t/;
    $h{ $f[0].$f[1].$f[2].$f[3] } = '  Ref='.$_;
  }
  return \%h;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
  print "Usage: $0 [options] nesoni.csq.1 nesoni.csq.2\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
