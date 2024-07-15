#!/usr/bin/env perl
use strict;
use warnings;
use Fatal;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fasta write_fastq assert);

my(@Options, $verbose, $num, $prop, $fasta);
setOptions();

die "please specify -n or -p" unless ($num or $prop);
die "bad -p $prop : must be between 0 and 1" if $prop and $prop < 0 or $prop > 1;

# open input file

my $Lfn = shift @ARGV or die "please specify input FASTQ file";
print STDERR "Reads: $Lfn\n";
-r $Lfn or die "can't read '$Lfn'";

my $cmd = $Lfn =~ m/gz$/i ? "zcat -f $Lfn |" : "< $Lfn";
print STDERR "Opening using '$cmd'\n";
open my $Lfh, $cmd;

# count number of sequences

print STDERR "Counting sequences: $Lfn ...\n";
my($nseq) = qx(zcat -f $Lfn | wc -l);
chomp $nseq;
die "$Lfn had $nseq lines - not a multiple of 4" unless $nseq % 4 == 0;
die "You chose $num seqs, but file only has $nseq" if $num > $nseq;
$nseq = $nseq/4;
$num = int($prop*$nseq) if $prop;
print STDERR "Found $nseq, choosing random $num.\n";

# set up selection bitmap

print STDERR "Building random selection bitmap...\n";
my @pick = map { 0 } (1 .. $nseq);
print STDERR "Looping... please press ENTER to seed the RNG \n";
my $picked=0;
while ($picked < $num) {
  my $i = int rand($nseq-1);
  next if $pick[$i];
  $pick[$i] = 1;
#  print STDERR "picked $i\n";
  $picked++;
}


# filter inputs via selection map

my $writer_func = $fasta ? \&write_fasta : \&write_fastq;

my $found=0;
my $index=0;
while ( not eof () ) { 
  my $L = read_fastq($Lfh);
  next unless $pick[$index++];
  $writer_func->(\*STDOUT, $L);
  $found++;
  print STDERR "\rSelected: $found";
  last if $found >= $num;
}
print STDERR "\nDone.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"verbose info"},
    {OPT=>"n|num=i",  VAR=>\$num, DEFAULT=>0, DESC=>'Choose N sequences from input'},
    {OPT=>"p|prop=f",  VAR=>\$prop, DEFAULT=>0, DESC=>'Choose this proportion of input sequences'},
    {OPT=>"f|fasta!",  VAR=>\$fasta, DEFAULT=>0, DESC=>'Output FASTA instead of FASTQ'},
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
  print "Usage: $0 [options] <file.fq[.gz]>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
