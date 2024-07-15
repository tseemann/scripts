#!/usr/bin/env perl
use strict;
use warnings;
use Fatal;
use List::Util qw(sum min max);
use Data::Dumper;
use IO::File;

my(@Options, $verbose);
setOptions();

my $maxL = 0;
my @sum;
my @cnt;
my @num;

for my $file (@ARGV) { 
  my $nseq=0;
  my $fh = open_maybe_compressed($file);
  while (my $id = $fh->getline) {
    print STDERR "\rReading: $nseq" if ++$nseq % 971 == 0;
    if ($id =~ m/^\@/) {
      $fh->getline; # skip dna
      $fh->getline; # skip id2
      my $qs = $fh->getline;
      # this is slow... maybe use unpack() or use C ?
      my @qv = split m//, $qs;
      @qv = map { ord } @qv;
      $maxL = max($maxL, scalar @qv);
      $num[$#qv]++; # count reads of this length
      for my $i (0 .. $#qv) {
        $sum[$i] += $qv[$i];
	$cnt[$i] += 1;
      }
    }
  }
  print STDERR "\rRead $nseq sequences from $file\n";
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
print STDERR "Guessed phred offset = $offset ($name_of{$offset})\n";

# print table
#print Dumper(\@sum, \@num) if $verbose;
for my $i (0 .. $maxL-1) {
  printf "%d\t%.1f\t%d\n", $i+1, $sum[$i]-$offset, ($num[$i] || 0);
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
      print STDERR "open_maybe_compressed: using $exe to read $fname\n";
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
