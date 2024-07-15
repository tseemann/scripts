#!/usr/bin/env perl
use strict;
use Data::Dumper;
use List::MoreUtils qw(each_array);

# "qseq"
#HWUSI-EAS-100R  0002    3       2       999     6334    0       1       ............................................................................    BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB    0

# "fastq"
#@ILLUMINA-48B6CB_0001:1:1:1057:3838#CGATGT/1
#CTTCTAGTAAAATGTATCTTCATATCGGCTTCCTCACTTTTCTTCCGTAACAGGATCAGAGAATTTCGGAATAAGAACTCCCTGCCCCACTGCGAAATTC

my(@Options, $verbose, $lane, $read, $destdir);
setOptions();

my $dir = shift @ARGV;
print STDERR "Folder: $dir\n";

my @lane = $lane ? ($lane) : (1..8);
print STDERR "Lanes: @lane\n";

$read >= 0 && $read <= 2 or die "--read 0|1|2 allowed only";
#$read = 3 if $read==2; # mux=2
my @read = $read ? ($read) : (1, 2);
print STDERR "Reads: @read\n";

for my $lane (@lane) {
  for my $read (@read) {
    print STDERR "Doing: Lane $lane, Read $read\n"; 
    my @seq = get_qseq_files($dir, $lane, $read==1 ? 1 : 3);
    printf STDERR "\tFound %d Read$read files.\n", scalar @seq;
    my @mux = get_qseq_files($dir, $lane, 2);
    printf STDERR "\tFound %d MuxTag files.\n", scalar @mux;
    die "different number of seq and mux files!" if @seq != @mux;
    my $dest = "$destdir/s_${lane}_${read}_sequence.txt";
    open my $dfh, '>', $dest;
    print "Writing: $dest\n";
    my $each = each_array(@seq, @mux);
    while ( my($sf,$mf) = $each->() ) {
      print STDERR "Processing: $sf + $mf\n";
      open my $sfh, '<', "$dir/$sf";
      open my $mfh, '<', "$dir/$mf";
      append_fastq($dfh, $sfh, $mfh);
    }
  }
}

#----------------------------------------------------------------------

sub append_fastq {
  my($dfh, $sfh, $mfh) = @_;
  while (my $seq = <$sfh>) {
    my $mux = <$mfh>;
    my @s = split m/\t/, $seq;
    next if substr($s[10],0,1) eq '0'; # failed purity test
    chop $s[8]; # last base is unphased or somewhat, so chop it
    $s[8] =~ s/\./N/g;
    my @m = split m/\t/, $mux;
    chop $m[8]; # last base in muxtag is for phasing (should be 'A')
    $m[8] =~ s/\./N/g;
    next if substr($m[10],0,1) eq '0'; # failed purity test
    my $read = $s[7] == 1 ? 1 : 2;
    print $dfh "\@$s[0]_$s[1]:", join(':',@s[2..5]),"#$m[8]/$read\n";
    print $dfh $s[8],"\n+\n",$s[9],"\n";
  }
}

#----------------------------------------------------------------------

sub get_qseq_files {
  my($dir, $lane, $read) = @_;
  opendir DIR, $dir;
  my @file = sort grep { m/s_${lane}_${read}_\d+_qseq.txt/ } readdir DIR;
  closedir DIR;
#  return @file[0..0];
  return @file;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"l|lane=i",  VAR=>\$lane, DEFAULT=>0, DESC=>"Just this lane"},
    {OPT=>"r|read=i",  VAR=>\$read, DEFAULT=>0, DESC=>"Just this read (1 or 2)"},
    {OPT=>"d|destdir=s",  VAR=>\$destdir, DEFAULT=>'.', DESC=>"Destination dir"},
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
  print "Usage: $0 [options] /bio/illumina/Runs/100506_HWUSI-EAS-100R_???????/Data/Intensities/BaseCalls/\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
