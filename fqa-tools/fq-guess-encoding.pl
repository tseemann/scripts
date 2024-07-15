#!/usr/bin/env perl
use strict;
use warnings;
use Fatal;
use List::Util qw(sum min max);
use Data::Dumper;
use IO::File;

my(@Options, $verbose, $numreads, $mode);
setOptions();

my $maxL = 0;
my @sum;
my @cnt;

FILE: 
for my $file (@ARGV) { 
  my $nseq=0;
  my $fh = open_maybe_compressed($file);
  while (my $id = $fh->getline) {
    last FILE if ++$nseq > $numreads;
    if ($id =~ m/^\@/) {
      $fh->getline; # skip dna
      $fh->getline; # skip id2
      my $qs = $fh->getline;
      my @qv = split m//, $qs;
      @qv = map { ord } @qv;
      $maxL = max($maxL, scalar @qv);
      for my $i (0 .. $#qv) {
        $sum[$i] += $qv[$i];
	$cnt[$i] += 1;
      }
    }
    else {
      die "File '$file' does not look like FASTQ";
    }
  }
}

# convert to average pseudo-Q values (need to subtract offset still)
for my $i (0 .. $maxL-1) {
  $sum[$i] ||= 0;
  $sum[$i] /= $cnt[$i];
}

# Guess ASCII encoding offset
# Sanger   = 33 # Solexa   = 59 # Illumina = 64
my %name_of = (33=>'Sanger', 59=>'Solexa', 64=>'Illumina');
my $minq = min(@sum);
my $offset = $minq < 59 ? 33 : ($minq < 65 ? 59 : 64);
print STDERR "Guessed phred offset = $offset ($name_of{$offset})\n" if $verbose;

my $stdout = '';
if ($mode eq 'offset') {
  $stdout  = $offset;
}
elsif ($mode eq 'bowtie2') {
  $stdout = "--phred-$offset";
}
elsif ($mode eq 'nesoni') {
  $stdout = "--qv-offset $offset";
}
elsif ($mode eq 'fqatools') {
  $stdout = "-o $offset";
}
elsif ($mode eq 'fastx') {
  $stdout = "-Q $offset";
}
else {
  die "unknown mode '$mode'";
}

print STDOUT "$stdout\n";

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
    {OPT=>"numreads=i",  VAR=>\$numreads, DEFAULT=>10_000, DESC=>"Sample this many reads"},
    {OPT=>"mode=s",  VAR=>\$mode, DEFAULT=>'offset', DESC=>"Mode: offset fastx bowtie2 nesoni fqatools"},
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
