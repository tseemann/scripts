#!/usr/bin/env perl
use strict;
use warnings;
use Fatal;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fasta write_fastq assert);
use Data::Dumper;
use List::Util qw(sum min max);

my(@Options, $verbose);
setOptions();

my %num;
my %runs;
my %runlen;
my %len;
my %goodlen;

while ( not eof() ) {
  print STDERR "\rReading: $." if $. % 19873 == 0;
  my $r = read_fastq( \*ARGV );
#  print $r->[1], "\n";
 
  my $L = length $r->[1];
  $len{$L}++;
 
  my $Ns = $r->[1] =~ tr/N//; # count Ns non-destructively
  $num{$Ns}++;
  
  my @Nrun = ($r->[1] =~ m/(N+)/gi);
#  print scalar(@Nrun), ": @Nrun\n";
  $runs{ scalar(@Nrun) } ++;

  my @good = ($r->[1] =~ m/[AGTC]+/gi);
  if (@good) {
    @good = map { length $_ } @good;
    @good = sort @good;
    $goodlen{ $good[-1] }++;
  }
  else {
    $goodlen{0}++;
  }

  for my $x (@Nrun) {
    $runlen{ length($x) }++;
  }
  
#  last if $. == 100;
}
my $nreads = sum(values %num);
print STDERR "\rProcessed $nreads reads.\n";

print STDERR Dumper(\%num) if $verbose;

my $L = max(keys %len);

print "VALUE,RAW READS THIS LONG,NON-N SUB-READS THIS LONG,READS WITH THIS MANY Ns,NO OF HOMO-N-POLYMERS PER READ,LENGTH OF HOMO-N-POLYMERS\n";
for my $i (0 .. $L) {
  print join(",",
    $i, 
    $len{$i} || 0,
    $goodlen{$i} || 0,
    $num{$i} || 0,
    $runs{$i} || 0,
    $runlen{$i} || 0,
  ), "\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose info"},
  );

  (@ARGV < 1) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] <reads.fastq>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
