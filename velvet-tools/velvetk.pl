#!/usr/bin/perl 
use strict;
use warnings;
use Fatal;
use IO::File;
use IO::Handle;
use List::Util qw(sum);
use Data::Dumper;

#    velvetk.pl             
#
#    Examine k-mer coverage of your reads for a given target genome size.
#    Assumes uniform coverage of your genome.
#
#    Copyright (C) 2012- Torsten Seemann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

my(@Options, $verbose, $size, $cov, $range, $genome, $best, $trimq, $qoffset);
setOptions();

die "Please only use one of --size or --genome" if $size and $genome;
if ($genome) {
  print STDERR "Estimating target genome size from '$genome'\n";
  $size=0;
  # doing this manually to avoid BioPerl dependency
  open my $in, '<', $genome;
  while (<$in>) {
    next if m/>/;
    chomp;
    $size += length $_;
  }
}

my %multiplier = (K=>1E3, M=>1E6, G=>1E9);
if ($size =~ m/(.*)([KMG])$/i) {
  $size = $1 * $multiplier{uc $2};
}

die "Target genome size '$size' isn't an integer" if $size !~ m/^[0-9E\-+]+$/i;
die "Target genome size is zero... can't proceed" if $size <=0;
print STDERR "Target genome size is $size bp\n";
print STDERR "Desire k-mer coverage of $cov\n" if $best;

my %num_with_len;
my $dna;
for my $file (@ARGV) { 
  my $nseq=0; 
  my $fh = open_maybe_compressed($file);
  while (my $line = $fh->getline) {
    print STDERR "\rReading: $nseq" if $nseq % 9871 == 0;
    if ($line =~ m/^>/) {
      # ok,fasta
      $dna = $fh->getline;
      chomp $dna;
    }
    elsif ($line =~ m/^\@/) {
      # ok, fastq
      $dna = $fh->getline;
      chomp $dna;
      $fh->getline; # skip id2
      my $qv = $fh->getline;
#      print $qv;
      chomp $qv;
      if (0 and $trimq > 0) {
        while (ord(substr $qv, -1, 1)-$qoffset < $trimq) {
          printf "'%d' %d %s %s\n", substr($qv,-1,1), ord(substr($qv,-1,1))-$qoffset, $dna, $qv;
          chop $dna;
          chop $qv;
        }
#        print "Trimmed: $dna\n";
      }
    }
    else {
      # dunno, skip and try next line!
      next;
    }
    # process the seq line
    $num_with_len{ length($dna) }++;
    $nseq++;
  }
  print STDERR "\rRead $nseq sequences from $file\n";
}

my @l = sort { $a <=> $b } keys %num_with_len;
die "no reads found" if not @l;

my $nreads = sum( values %num_with_len );
print STDERR "Considered $nreads reads with lengths $l[0]..$l[-1] bp\n";
print Dumper(\%num_with_len) if $verbose;

my $bestkc = 0;
my $bestk = 0;

print "K\t#Kmers\tKmer-Cov\n" if not $best;
for (my $k=1; $k <= $l[-1]; $k+=2) {
  my $n = num_kmers($k, \%num_with_len);
  my $kc = $n/$size;
  if (abs($kc-$cov) <= abs($bestkc-$cov)) {
    $bestk = $k;
    $bestkc = $kc;
  }
#  next if $kc < $cov-$range-0.5 or $kc > $cov+$range+0.5;
  printf "$k\t$n\t%.1f\n", $kc if not $best;
}
print "$bestk\n" if $best;

#----------------------------------------------------------------------

sub num_kmers {
  my($k, $num_with_len) = @_;
  my $n=0;
  for my $L (keys %{$num_with_len}) {
    next if $L < $k;
    $n += ($L-$k+1) * $num_with_len->{$L};
  }
  return $n;
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
      $io->ungetc(ord $c);
      print STDERR "Extracting reads via '$exe $fname |'\n";
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
    {OPT=>"size=s",  VAR=>\$size, DEFAULT=>0, DESC=>"Target genome size in bp (can use k/M/G suffix)"},
    {OPT=>"genome=s",  VAR=>\$genome, DEFAULT=>'', DESC=>"Target genome (FASTA)"},
    {OPT=>"cov=i",  VAR=>\$cov, DEFAULT=>25, DESC=>"Target k-mer coverage"},
#    {OPT=>"range=i",  VAR=>\$range, DEFAULT=>10, DESC=>"Print target +/- this range"},
    {OPT=>"best!",  VAR=>\$best, DEFAULT=>0, DESC=>"Just print best k-mer to stdout"},
#    {OPT=>"trimq=i",  VAR=>\$trimq, DEFAULT=>10, DESC=>"Trim FASTQ reads from 3' end < this Q"},
#    {OPT=>"qoffset=i",  VAR=>\$qoffset, DEFAULT=>33, DESC=>"FAST encoding: Sanger=33 Illumina<2012=64"},
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
  print STDERR "Synopsis:\n  List suitable k-mer sizes given target genome and reads\n";
  print STDERR "Author:\n  Torsten Seemann <torsten\@seemann.id.au>\n";
  print STDERR "Usage:\n  $exe [--size X | --genome F] [options] <readfile[.gz][.bz2] ...>\n";
  print STDERR "Options:\n";
  foreach (@Options) {
    printf STDERR "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
