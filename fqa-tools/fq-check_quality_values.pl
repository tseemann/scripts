#!/usr/bin/env perl
use strict;
use warnings;
use Fatal;
use List::Util qw(sum min max);
use Data::Dumper;
use IO::File;

my(@Options, $verbose, $numreads, $qvoffset);
setOptions();

my @qval;
my $minq=256;
my $maxq=-1;
my $nseq=0;

FILE: 
for my $file (@ARGV) { 
  my $fh = open_maybe_compressed($file);
  while (my $id = $fh->getline) {
    last FILE if ++$nseq > $numreads;
    print STDERR "\rProcessing: $nseq" if $nseq % 10_000 == 0;
    if ($id =~ m/^\@/) {
      $fh->getline; # skip dna
      $fh->getline; # skip id2
      my $qs = $fh->getline;
      chomp $qs;
      my @qv = split m//, $qs;
      @qv = map { ord } @qv;
      $qval[$_]++ for @qv;  # do with hash slice?
      $minq = min($minq, min(@qv));
      $maxq = max($maxq, max(@qv));
    }
    else {
      die "File '$file' does not look like FASTQ";
    }
  }
}
print STDERR "Processed $nseq sequences.\n";

print "Symbol\tASCII\tQ+$qvoffset\tFrequency\n";
for my $Q ($minq .. $maxq) {
  my $sym = $Q < 32 ? '<np>' : chr($Q);
  print join("\t", $sym, $Q, $Q-$qvoffset, $qval[$Q]||0), "\n";
}

#----------------------------------------------------------------------

sub open_maybe_compressed {
  my($fname) = @_;
  my @try = ("pbzip2 -dc", "pigz -dc", "bzip2 -dc", "gzip -dc", "cat");
  my $fh;
  for my $exe (@try) {
    print STDERR "Trying: $exe $fname\n" if $verbose;
    my $io = IO::File->new("$exe \Q$fname\E 2> /dev/null |");
    my $c = $io->getc;
    if (defined $c) {
      $io->ungetc(ord $c); # is ord() needed here?
      print STDERR "open_maybe_compressed: using $exe to read $fname\n" if $verbose;
      return $io; 
    }
  }
  die "could not open $fname";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"numreads=i",  VAR=>\$numreads, DEFAULT=>1E20, DESC=>"Sample this many reads"},
    {OPT=>"qvoffset=i",  VAR=>\$qvoffset, DEFAULT=>33, DESC=>"Quality value offset"},
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
  my $exe = $0;
  $exe =~ s{^.*/}{};
  print STDERR "Usage: $exe <readfile[.gz][.bz2] ...>\n";
  foreach (@Options) {
    printf STDERR "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
