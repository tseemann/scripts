#!/usr/bin/env perl
use strict;
use Data::Dumper;
use IO::File;
use List::Util qw(min max);
use List::MoreUtils qw(uniq);
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fastq write_fasta assert);

my(@Options, $verbose, $prefix, $muxes, $fasta, $samplesize, $badmux);
setOptions();

assert(@ARGV >= 2, "please provide two .fq files as input");
assert(-r $ARGV[0], "can not open left reads: $ARGV[0]");
assert(-r $ARGV[1], "can not open right reads: $ARGV[1]");
#assert($muxes > 1, "number of muxed samples must be more than one");

$ARGV[0] =~ m/(\d)/; # first numeral in filename... usually works
my $lane = $1 || '0';
print STDERR "Detected lane source as: $lane\n";

open my $L, '<', $ARGV[0];
open my $R, '<', $ARGV[1];

printf STDERR "Multiplex: %s\n", $muxes || 'AUTO-DETECT';

my @barcode = get_barcodes($ARGV[0]);
my $nbc = scalar @barcode;
print STDERR "Detected $nbc barcodes/muxtags\n";
if ($nbc == 1) {
  print STDERR "Only $nbc barcode - no need to demux!\n";
  exit;
}
print STDERR "Using barcodes: @barcode\n";
my %muxid_of = map { ($barcode[$_] => $_+1) } (0..$#barcode);
#exit;

for my $bc (@barcode) {
  my @sim = close_matches($bc);
  
  if ($badmux) {
    my @more;
    for my $sim (@sim) {
      push @more, close_matches($sim);
    }
    @more = uniq(@more);
    @sim = @more;
  }
  
  print STDERR "Mapping: $bc <= @sim\n" if $verbose;
  map { $muxid_of{$_} = $muxid_of{$bc} } @sim;
}
#print Dumper(\%muxid_of);
#exit;

my %fh_of;
for my $bc (@barcode) {
  for my $dir (1, 2) {
    my $name = $prefix;
    $name =~ s/\%l/$lane/g;
    $name =~ s/\%b/$bc/g;
    $name =~ s/\%s/$muxid_of{$bc}/g;
    $name =~ s/\%d/$dir/g;
    print STDERR "Opening for writing: $name\n";
    $fh_of{$muxid_of{$bc}}{$dir} = IO::File->new(">$name");
  }    
}

print STDERR "Demuxing: $ARGV[0] $ARGV[1]\n";
my $nread=0;
my $nwrote=0;

my $writer_func = $fasta ? \&write_fasta : \&write_fastq;

while (not eof $L) {
  my $l = read_fastq($L);
  my $r = read_fastq($R);
  $nread++;

  next unless $l->[0] =~ m/\#([AGTC]+)\/\d$/;
  my $bc = $1;
  my $mid = $muxid_of{$bc};
  next unless $mid;
                       
  $nwrote++;
  
  $writer_func->($fh_of{$mid}{1}, $l);
  $writer_func->($fh_of{$mid}{2}, $r);
}

printf STDERR "Parsed %d pairs. Discarded %d. Wrote %d (%.2f%%)\n", 
  $nread, 
  $nread-$nwrote, 
  $nwrote, ($nwrote*100/$nread);

#----------------------------------------------------------------------

sub get_barcodes {
  my($lfn) = @_;
  my %freq;
  my $remaining = $samplesize;
  print STDERR "Guessing barcodes from $remaining reads in $lfn ...\n";
  open my $lfh, '<', $lfn;
  while (not eof $lfh) {
    my $r = read_fastq($lfh);
    $r->[0] =~ m/\#([NAGTC]+)\/\d$/;
    $1 or die "No mux tag or bar code found in read ID: ".$r->[0];
    next if $1 =~ m/N/;
    $freq{$1}++;
    last if --$remaining <= 0;
  }
  close $lfh;
  my $M = guess_mux(\%freq);
  my @seen = sort { $freq{$b} <=> $freq{$a} } keys %freq;
  return @seen[0 .. $M-1];
}

#----------------------------------------------------------------------

sub guess_mux {
  my($freq) = @_;
  # if --muxes was chosen
  return $muxes if $muxes;
  # else auto-detect;
  print STDERR Dumper($freq) if $verbose;
  my $max = max(values %$freq);
  return scalar( grep { $_ > $max/2 } values %$freq);
}

#----------------------------------------------------------------------
# generate all barcodes that are 1 SNP away from the reference

sub close_matches {
  my($ref) = @_;
  my @list;
  for my $i (0 .. length($ref)-1) {
    my $base = substr($ref,$i,1);
    for my $sub (grep { $_ ne $base } qw(A G T C N)) {
      my $new = $ref;
      substr($new,$i,1) = $sub;
      push @list, $new;
    }
  }
  return @list;
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"p|prefix=s", VAR=>\$prefix, DEFAULT=>'s_%l_%d_%b.txt', DESC=>"Output naming: %l=lane %b=barcode %s=muxID %d=pairID"},
    {OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
    {OPT=>"m|muxes=i", VAR=>\$muxes, DEFAULT=>0, DESC=>"Number of multiplexed samples (0=autodetect)"},
    {OPT=>"s|samplesize=i", VAR=>\$samplesize, DEFAULT=>500000, DESC=>"Number of reads to guess mux tags from"},
    {OPT=>"b|badmux!", VAR=>\$badmux, DEFAULT=>0, DESC=>"Allow 2 errors in barcode (lots of Ns...)"},

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
  print "Usage: $0 [--muxes N] left.fq right.fq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
