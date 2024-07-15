#!/usr/bin/env perl
use strict;

my(@Options, $verbose, $slow, $exact, $noheader, $name);
setOptions();

@ARGV = map { m/\.(gz|Z)$/ ? "gzip -dc < $_ |" : $_ } @ARGV;

my $count=0;
my $letters=0;
my $minlen=1E20;
my $maxlen=0;

sub fast {
while (<>) { # ID
  $_ = <>;   # sequence
  chomp;
  my $L = length;
  $count++;
  $letters += $L;
  $minlen = $L if $L < $minlen;
  $maxlen = $L if $L > $maxlen;
  $_ = <>; # ID again
  $_ = <>; # quality
}
}

sub slow {
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq assert);
use List::Util qw(max min);
while ( not eof() ) {
  my $s = read_fastq( \*ARGV );
  my $L = length($s->[1]);  # 3-tuple [ ID,seq,qual ]
  $count++;
  $letters += $L;
  $minlen = min($minlen, $L);
  $maxlen = max($maxlen, $L);
}
}

$slow ? slow() : fast();

my @COL = qw(reads bp minlen maxlen);
unshift @COL, 'name' if $name;

print join ("\t", @COL),"\n" unless $noheader;

my @res = ($count, $letters, $minlen, $maxlen);
@res = map { downsize($_) } @res unless $exact;
unshift @res, $name if $name;
print join ("\t", @res),"\n";

sub downsize {
  my($n) = @_;
  if (length($n) > 9) {
    return int($n/1E9).'G';
  }
  elsif (length($n) > 6) {
    return int($n/1E6).'M';
  }
  elsif (length($n) > 3) {
    return int($n/1E3).'k';
  }
  else {
    return $n;
  }
  
}      

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"slow!", VAR=>\$slow, DEFAULT=>0, DESC=>"Use slow function calls"},
    {OPT=>"exact!", VAR=>\$exact, DEFAULT=>0, DESC=>"Output exact numbers, not k/M/G approximations"},
    {OPT=>"noheader!", VAR=>\$noheader, DEFAULT=>0, DESC=>"Don't output header"},
    {OPT=>"name=s", VAR=>\$name, DEFAULT=>'', DESC=>"If set, add column with this name"},
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
  print "Usage: $0 [--each] file1.fq [ file2.fq.gz file3.fq.Z ... ]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
