#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my(@Options, $verbose, $folder);
setOptions();

my $ref_fn = shift @ARGV;
die "can't read '$ref_fn' reference gbk/fasta" unless -r $ref_fn;
my($line) = qx(head -1 $ref_fn);
die "reference doesn't look like .gbk or .fasta" unless $line =~ m/(^>|^LOCUS)/;
print STDERR "Reference: $ref_fn\n";

my @dir = @ARGV;
die "please specify some Dropbox folders." unless @dir;
print STDERR "Dropbox folders: ", 0+@dir, "\n";
my $count=0;
my @sample;
for my $dir (@dir) {
  -d $dir or die "'$dir' is not a valid folder";
  my($s) = <$dir/analysis/s_*.fa>;
  $s && -r $s or die "Can't fine s_*.fa in $dir/analysis/";
  my $size = -s $s;
  $size or die "$s is an empty file.";
  $size = int $size/1024/1024;
  $count++;
  $dir =~ m{([^/]+)/?$};
  die "can't figure out sample name from '$dir'" unless $dir;
  my $name = $1;
  print STDERR "[$count] $name | $s | $size MB\n";
  push @sample, [ $name, $s ];
}
print STDERR "Samples: ", 0+@sample, "\n";
print "mkdir -p '$folder'\n\n";

for my $sample (@sample) {
  my($name, $reads) = @$sample;
  my $ndir = "$folder/$name";
  print "\n\n# $name\n\n";

  my $cmd = "pypy-nesoni samshrimp: $ndir".
            " --sam-unaligned no".
	    " $ref_fn".
	    " reads: $reads".
	    " shrimp-options: -h 68%".
	    "";
  print "if [ ! -r '$ndir/alignments_filtered_sorted.bam.bai' ]; then\n";
  print "\t$cmd\n";
  print "fi\n"; 
  print "echo '$name shrimp done.'\n";
  
  $cmd = "pypy-nesoni samconsensus: $ndir".
         " --trim 3 --cutoff 0.9 --majority 0.75 --ambiguity-codes no".
	 "";
  print "if [ ! -r '$ndir/consensus.fa' ]; then\n";
  print "\t$cmd\n";
  print "fi\n"; 
  print "echo '$name consensus done.'\n";
  print "echo -n 'Number of SNPs found: '\n";
  print "tail --lines=+2 $ndir/report.txt | wc -l\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"f|folder=s",  VAR=>\$folder, DEFAULT=>'.', DESC=>"Where to put results"},
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
  print "Usage: $0 [options] <ref.gbk|ref.fna> <DropboxDir1> <DropboxDir2> ...\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
