#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fasta write_fasta assert);

my(@Options, $verbose, $newid, $start, $nons, $trunc, $nopoly, $chop, $minlen);
setOptions();

assert(@ARGV >= 2, "please provide two .fasta files as input");
assert(-r $ARGV[0], "can not open left reads: $ARGV[0]");
assert(-r $ARGV[1], "can not open right reads: $ARGV[1]");

open my $L, '<', $ARGV[0];
open my $R, '<', $ARGV[1];

if ($verbose) {
  print STDERR "Interleaving $ARGV[0] + $ARGV[1]\n";
  print STDERR "Removing pairs with 'N's.\n" if $nons;
  print STDERR "Truncating reads to $trunc bp.\n" if $trunc;
  print STDERR "Chopping last base off reads\n" if $chop;
  print STDERR "Please be patient...\n";
}

my $nread=0;
my $nwrote=0;

while (not eof $L) {
  my $l = read_fasta($L);
  my $r = read_fasta($R);
  $nread++;
  print STDERR "\rProcessed: $nread" if $nread % 9751 == 0;

  # chop
  if ($chop) {
    chop $l->[1];
    chop $r->[1];
  }
  
  # start
  if ($start > 0) {
    $l->[1] = substr($l->[1], $start);
    $r->[1] = substr($r->[1], $start);
  }

  # truncate 
  if ($trunc > 0) {
    $l->[1] = substr($l->[1], 0, $trunc);
    $r->[1] = substr($r->[1], 0, $trunc);
  }

  # then skip read if either left or right has N in it
  next if $nons and index($l->[1].$r->[1], 'N') >= $[ ;
  
  # don't allow monopolymers through ?
  next if $nopoly and ( $l->[1] =~ m/^(.)\1*$/i or $r->[1] =~ m/^(.)\1*$/i );
                       
  next if length($l->[1]) < $minlen or length($r->[1]) < $minlen;
		       
  $nwrote++;

  if ($newid) {
    $l->[0] = "$nwrote/1";
    $r->[0] = "$nwrote/2";
  }

  write_fasta(\*STDOUT, $l);
  write_fasta(\*STDOUT, $r);
}


printf STDERR "\rParsed %d pairs. Discarded %d. Wrote %d (%.2f%%)\n", 
  $nread, 
  $nread-$nwrote, 
  $nwrote, ($nwrote*100/$nread);


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"i|newid!", VAR=>\$newid, DEFAULT=>0, DESC=>"Generate new, compact IDs"},
    {OPT=>"n|nons!", VAR=>\$nons, DEFAULT=>0, DESC=>"Skip pairs with Ns"},
    {OPT=>"p|nopoly!", VAR=>\$nopoly, DEFAULT=>0, DESC=>"Skip monopolymer reads"},
    {OPT=>"chop!", VAR=>\$chop, DEFAULT=>0, DESC=>"Chop last base (3' end)"},
    {OPT=>"t|trunc=i", VAR=>\$trunc, DEFAULT=>0, DESC=>"Truncate reads to this length (0 = no truncation)"},
    {OPT=>"s|start=i", VAR=>\$start, DEFAULT=>0, DESC=>"Trim read to this start coordinate (0 = ignore)"},
    {OPT=>"m|minlen!", VAR=>\$minlen, DEFAULT=>21, DESC=>"Minimum length to allow through"},
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
  print "Usage: $0 [--fasta] left.fa right.fa > interleaved.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
