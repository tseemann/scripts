#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fastq trim_fastq assert);

my(@Options, $verbose, $fasta, $minq, $minlen, $suffix, $prefix);
setOptions();

assert($minq >= -5 and $minq <= 40, "--minq must be between -5 and 40");
assert(@ARGV >= 2, "please provide two .fq files as input");
assert(-r $ARGV[0], "can not open left reads: $ARGV[0]");
assert(-r $ARGV[1], "can not open right reads: $ARGV[1]");

open my $L, '<', $ARGV[0];
open my $R, '<', $ARGV[1];
open my $Lout, '>', "$prefix$ARGV[0]$suffix";
open my $Rout, '>', "$prefix$ARGV[1]$suffix";

#print STDERR "# trimming bases with Q < $minq\n";
#print STDERR "# rejecting trimmed reads containing an N\n" if not $allown;
print STDERR "Trimming $ARGV[0] $ARGV[1] pairs to Q$minq. Please be patient...\n";

my($nread,$ntrim,$nambig,$nshort,$npass, $oldbp,$newbp) = (0)x7;

while (not eof $L) {
  my $l = read_fastq($L);
  my $r = read_fastq($R);
  $nread++;
  $oldbp += length($l->[1]) + length($r->[1]);
  
  $l = trim_fastq($l, $minq);
  next unless defined $l;
  $r = trim_fastq($r, $minq);
  next unless defined $r;
  $ntrim++;
  
  next if index($l->[1],'N') > 0 or index($r->[1],'N') > 0;
  $nambig++;
  
  next if length($l->[1]) < $minlen or length($r->[1]) < $minlen;
  $nshort++;
  
  write_fastq($Lout, $l);
  write_fastq($Rout, $r);
  $npass++;
  $newbp += length($l->[1]) + length($r->[1]);
}

if ($nread) {
  printf STDERR "%d reads, trimmed(Q$minq) %d, unambig %d, size(${minlen}bp) %d, passed %d (%.2f%%)\n",
    $nread, $ntrim, $nambig, $nshort, $npass, $npass*100/$nread;
  printf STDERR "read %d bp, wrote %d bp (%.2f%%)\n",
    $oldbp, $newbp, $newbp*100/$oldbp;
}
else {
  print STDERR "No reads found.\n";
}    


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
#    {OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
    {OPT=>"minq=i", VAR=>\$minq, DEFAULT=>15, DESC=>"Minimum quality: -5 to 40"},
    {OPT=>"minlen=i", VAR=>\$minlen, DEFAULT=>31, DESC=>"Minimum read length in output"},
    {OPT=>"prefix=s", VAR=>\$prefix, DEFAULT=>'', DESC=>"Prefix for output files"},
    {OPT=>"suffix=s", VAR=>\$suffix, DEFAULT=>'.trim', DESC=>"Suffix for output files"},
#    {OPT=>"s|singletons=s", VAR=>\$singletons, DEFAULT=>'.singeltons', DESC=>"Suffix for singletons"},
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
  print "Usage: $0 [options] left.fq right.fq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
