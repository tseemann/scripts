#!/usr/bin/env perl
use strict;
use warnings;
use Fatal;
use File::Temp qw(tempdir);
use File::Path qw(remove_tree);

my(@Options, $quiet, $nodesc, $fasta, $tempdir, $rev, $min, $max);
setOptions();

die "ERROR: --max ($max) is less than --min ($min)" if $max < $min;
my $verbose = ! $quiet;

my $dir;      # temp dir
my @line;     # a fastq entry
my @fh;       # array of output filehandles, one per length (radix sort)
my $nread=0;  # num read
my $nwrote=0; # num wrote

$SIG{INT} = sub { 
    remove_tree($dir) if defined $dir; 
};

while ( ! eof () ) {
  $line[$_] = scalar(<>) for (0..3);
  $nread++;
  print STDERR "Passed $nwrote/$nread\n" if $verbose and ($nread % 100000 == 0);
  
  my $len = length($line[1]) - 1;
  next if $len < $min;
  next if $len > $max;
  
  if ( not defined $fh[$len] ) {
    if (not defined $dir) {
      # only create temp dir once we start writing
      $dir = tempdir(DIR=>$tempdir, CLEANUP=>1);
      print STDERR "Using temp folder: $dir\n" if $verbose;      
    }
    open $fh[$len], '>', "$dir/$len";
    print STDERR "Creating file: $dir/$len\n";
  }
  
  if ($nodesc) {
    $line[0] =~ s/\s.*$//;
  }
  
  if ($fasta) {
    substr $line[0], 0, 1, '>'; # replace @ with >
    print {$fh[$len]} $line[0], $line[1];
  }
  else {
    print {$fh[$len]} @line;
  }
  $nwrote++;
  
}
print STDERR "Writing out sorted files...\n" if $verbose;

my @order = (0 .. $#fh);
@order = reverse @order if $rev;

for my $len (@order) {
  next unless defined $fh[$len];
  system("cat", "$dir/$len");
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"quiet!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Quiet (no output)"},
    {OPT=>"tempdir=s",  VAR=>\$tempdir, DEFAULT=>'.', DESC=>"Temp directory"},
#    {OPT=>"rename=s",  VAR=>\$rename, DEFAULT=>'', DESC=>"Rename reads with printf format string eg. read%09d"},
    {OPT=>"v|reverse!",  VAR=>\$rev, DEFAULT=>0, DESC=>"Reverse sort (long to short)"},
    {OPT=>"nodesc!",  VAR=>\$nodesc, DEFAULT=>0, DESC=>"Remove anything after first space in read name"},
    {OPT=>"fasta!",  VAR=>\$fasta, DEFAULT=>0, DESC=>"Output FASTA instead of FASTQ)"},
    {OPT=>"min=i",  VAR=>\$min, DEFAULT=>1, DESC=>"Minimum read length to keep"},
    {OPT=>"max=i",  VAR=>\$max, DEFAULT=>1E9, DESC=>"Maximum read length to keep"},
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
