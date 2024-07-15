#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my(@Options, $verbose, $project, $centre, $runfile, $srcdir, $outdir, $final);
setOptions();

my $SEP = ',';
#my $SEP = "\t";

my %run;
my $count=0;

-r $runfile or die "can't open --runfile $runfile";

mkdir $outdir or die "can't create folder: $outdir";
open SHEET2, ">$outdir/sheet2.csv";
open SHEET3, ">$outdir/sheet3.csv";
open SHEET4, ">$outdir/sheet4.csv";

open RUN, $runfile;
while (my $name = <RUN>) {
  chomp $name;
  print STDERR "Read: $name\n";
  next if $name =~ m/^\s*#/;
  my $dir = "$srcdir/$name";
  print STDERR "Folder: $dir\n";
  die "no such folder: $dir" unless -d $dir;
  $run{FOLDER} = $dir;
  $run{SAMPLEID} = sprintf "SAMNxxxx%04d", ++$count; 

  print SHEET2 join($SEP,
    $project,
    $run{SAMPLEID},
    $name,
  ), "\n";
  
  my @file = glob("$dir/*");
#  print STDERR "FILES: @file\n";

  
  my($platform,$machine,$layout,$len);
  my($sff_md5, $R1_md5, $R2_md5);
  
  my($sff) = grep { m/sff$/ } @file;
  my($pgm) = grep { m/\.fq.gz$/ } @file;  #HACK FOR MISSING SFF FILES
  
  if ($sff) {
    ($platform,$machine,$layout) = ('ION_TORRENT', 'Ion Torrent PGM', 'Single');
    print STDERR "SFF=$sff\n";
    system("ln -s \Q$sff\E $outdir/$name.sff");
    $sff_md5 = checksum($sff);
    $len = 200;
  }
  elsif ($pgm) {
    ($platform,$machine,$layout) = ('ION_TORRENT', 'Ion Torrent PGM', 'Single');
    print STDERR "PGM_FASTQ=$pgm\n";
    system("ln -s \Q$pgm\E $outdir/$name.fq.gz");
    $sff_md5 = checksum($pgm);
    $len = 200;
  }
  else {
    ($platform,$machine,$layout,$len) = ('ILLUMINA', 'Illumina HiSeq 2000', 'Paired-end');
    my($read1) = grep { m/_R1.*gz$|_1_sequence.txt.gz$/ } @file;
    $read1 or die "$name: can't find R1";
    my($read2) = grep { m/_R2.*gz$|_2_sequence.txt.gz$/ } @file;
    $read2 or die "$name: can't find R2";
    print STDERR "R1=$read1\nR2=$read2\n";
    system("ln -s \Q$read1\E $outdir/${name}_R1.fastq.gz");
    system("ln -s \Q$read2\E $outdir/${name}_R2.fastq.gz");
    $R1_md5 = checksum($read1);
    $R2_md5 = checksum($read2);
    $len = fastq_len($read1);
  }
  
  my $strain = $name;
  $strain =~ s/^.*?_//;
  
  print SHEET3 join($SEP, 
    "exp-$name",
    "Genome sequencing of Enterococcus faecium strain $strain",
    $centre,
    "lib-$name",
    $project,
    $run{SAMPLEID},
    'WGS',
    'GENOMIC',
    'size fractionation',
    $layout,
    $platform,
    $machine,
    'PLEASE FILL THIS IN',
    $len,
  ), "\n";

  
  print SHEET4 join($SEP, 
    "run-$name",
    "exp-$name",
    $centre,
    ($sff ? ('sff', "${name}.sff", $sff_md5)
          : $pgm ? ('fastq', "${name}.fq.gz", $sff_md5) :
                  ('fastq', "${name}_R1.fastq.gz", $R1_md5, 'fastq', "${name}_R2.fastq.gz", $R2_md5) ),
  ), "\n";

}

sub fastq_len {
  my($fname) = shift;
  print STDERR "Measuring FASTQ read length on: $fname\n";
  my($len) = qx(zcat \Q$fname\E | head -2 | tail -1);
  $len = length($len)-1;
  print STDERR "len: $len\n";
  return $len;
}

sub checksum {
  my($fname) = shift;
  print STDERR "Running md5sum on: $fname\n";
  my($res) = $final ? qx(md5sum \Q$fname\E) : ('fakemd5sumhere');
  chomp $res;
  $res =~ s/ .*$//; #chop off filename
  print STDERR "md5: $res\n";
  return $res;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"project=s",  VAR=>\$project, DEFAULT=>'PRJNA205886', DESC=>"Bioproject ID"},
    {OPT=>"centre=s",  VAR=>\$centre, DEFAULT=>'UOM', DESC=>"Submitting centre"},
    {OPT=>"runfile=s",  VAR=>\$runfile, DEFAULT=>'samples.txt', DESC=>"Text file of Dropbox folder names"},
    {OPT=>"srcdir=s",  VAR=>\$srcdir, DEFAULT=>'/bio/illumina/Dropbox/tim.stinear/E.faecium', DESC=>"Dropbox root folder"},
    {OPT=>"outdir=s",  VAR=>\$outdir, DEFAULT=>'sra-upload', DESC=>"Output folder"},
    {OPT=>"final!",  VAR=>\$final, DEFAULT=>0, DESC=>"Really do the md5 checksums! (slow)"},
  );

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
