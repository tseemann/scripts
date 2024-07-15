#!/usr/bin/env perl
use warnings;
use strict;
use Time::Piece;
use File::Copy; # cross-device safe file move instead of rename()

my(@Options, $verbose, $reference, $indels, $length, $single, $prefix, $coverage, $ambigs, $gzip);
setOptions();

msg("Checking 'mason' is installed");
my($mason) = qx(mason 2>&1);
$mason =~ m/read simulator/i or err("Please install the SeqAn 'mason' read simulator");

$reference or err("Please supply --reference genome.fasta");
-r $reference or err("Can't read reference: $reference");

#$prefix or err("Please provide an output filename prefix");
$prefix ||= $reference;
msg("Using '$prefix' prefix for output files");

my $logfile = "$prefix.log";
msg("Logging mason output to $logfile");

my $genomesize = -s $reference; # close enough, size in bytes
my $wantbases = int( $genomesize * $coverage );
my $numreads = int( ($single ? 1 : 0.5) * $wantbases / $length );
msg("$genomesize bp @ $coverage-fold coverage with $length bp reads = $numreads reads");

my $idp = $indels ? 0.00005 : 0;
my $pair = $single ? "" : "-mp -ll 500 -le 100";
my $psnp = 0.001;
my $enns = $ambigs ? "" : "-hnN -nN";

my $cmd = "mason illumina -sq -n $length -o mason.fasta -N $numreads -rn 2"
         ." -hs $psnp -hi $idp -hm 1 -hM 1"
         ." -pmm $psnp -pi $idp -pd $idp"
         ." -qmb 37 -qsdb 5  -qme 30 -qsde 10 -mmqmb 38 -mmqme 10"
         ." -seed $$ $enns $pair $reference 2> $logfile";

msg("Simulating reads, please be VERY patient...");
system($cmd)==0 or err("Error running mason - check $logfile for error");
msg("Success!");

msg("Removing temp files");
unlink "mason.fasta.sam";

msg("Cleaning up output files");
if ($single) {
  move("mason.fasta", "$prefix.fastq");
  msg("Reads are in: $prefix.fastq");
}
else {
  for my $i (1,2) {
    my $dest = $prefix."_R$i.fastq";
    move("mason_$i.fasta", $dest);
    msg("Read $i is in: $dest");
  }
}



msg("Done.");

#----------------------------------------------------------------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print STDERR $line;
}
        
#----------------------------------------------------------------------
        
sub err {
  msg(@_);
  exit(-1);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"prefix=s",  VAR=>\$prefix, DEFAULT=>'', DESC=>"Output filename prefix (.fastq or _R1.fastq suffix) [auto]"},
    {OPT=>"reference=s",  VAR=>\$reference, DEFAULT=>'', DESC=>"Reference genome (FASTA)"},
    {OPT=>"coverage=f",  VAR=>\$coverage, DEFAULT=>30, DESC=>"Target genome coverage / read depth"},
    {OPT=>"length=i",  VAR=>\$length, DEFAULT=>150, DESC=>"Read length"},
    {OPT=>"single!",  VAR=>\$single, DEFAULT=>0, DESC=>"Single end reads instead of paired-end"},
    {OPT=>"indels!",  VAR=>\$indels, DEFAULT=>0, DESC=>"Allow indels, not just subs"},
    {OPT=>"ambigs!",  VAR=>\$ambigs, DEFAULT=>0, DESC=>"Allow Ns in reads"},
#    {OPT=>"gzip!",  VAR=>\$gzip, DEFAULT=>0, DESC=>"Gzip output read files"},
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

__DATA__

mason - A Read Simulator
========================

SYNOPSIS
    mason illumina [OPTIONS] SEQUENCE

DESCRIPTION
    Use 'random' for the SEQUENCE file
    name to generate it randomly.

    -h, --help
          Displays this help message.
    --version
          Display version information

  Main Options:
    -aNg, --allow-N-from-genome
          Allow N from genome.
          Default: false.
    -s, --seed INT
          The seed for Rng. Default:
          0.
    -N, --num-reads NUM
          Number of reads (or mate
          pairs) to simulate. Default:
          1000.
    -sn, --source-length NUM
          Length of random source
          sequence. Default: 1000000.
    -spA, --source-probability-A NUM
          Propabilibty for A in
          randomly generated sequence.
          Default: 0.25.
    -spC, --source-probability-C NUM
          Propabilibty for C in
          randomly generated sequence.
          Default: 0.25.
    -spG, --source-probability-G NUM
          Propabilibty for G in
          randomly generated sequence.
          Default: 0.25.
    -snN, --source-no-N
          If set then no Ns are
          generated in the random
          source sequence.
    -f, --forward-only
          Simulate from forward strand
          only.
    -r, --reverse-only
          Simulate from reverse strand
          only.
    -o, --output-file STR
          Write results to PARAM.fasta
          file instead of
          SEQUENCE.reads.fasta.
    -sq, --simulate-qualities
          Simulate qualities, generate
          FASTQ instead of FASTA.
    -i, --include-read-information
          Include additional read
          information in reads file.
    -v, --verbose
          Verbosity mode.
    -vv, --very-verbose
          High verbosity mode, implies
          verbosity mode.
    -scf, --sample-counts-file STR
          Path to TSV file that maps
          contig ids to read counts.
    -vcf, --vcf-allele-file STR
          Path to VCF file to
          synthesize the haplotype
          from.

  Mate-Pair Options:
    -ll, --library-length-mean NUM
          Mate-pair mean library
          length. Default: 1000.
    -le, --library-length-error NUM
          Mate-pair library tolerance.
          Default: 100.
    -lu, --library-length-uniform
          Mate-pair library length is
          uniform.
    -mp, --mate-pairs
          Enable mate pair simulation.
    -rn, --read-naming NUM
          Read naming scheme in
          FASTQ/FASTA files. See
          section on read naming
          below. In range [0..2].
          Default: 0.
    -rnp, --read-name-prefix STR
          Read name prefix. Default is
          output file name.

  Haplotype Options:
    -hn, --num-haplotypes NUM
          Number of haplotypes to
          simulate. Default: 1.
    -hs, --haplotype-snp-rate NUM
          Haplotype SNP rate. Default:
          0.001.
    -hi, --haplotype-indel-rate NUM
          Haplotype indel rate.
          Default: 0.001.
    -hm, --haplotype-indel-range-min NUM
          Haplotype indel size min.
          Default: 1.
    -hM, --haplotype-indel-range-max NUM
          Haplotype indel size max.
          Default: 6.
    -hnN, --haplotype-no-N
          Do not allow Ns to be
          substituted or inserted into
          N. Default: Is allowed.

  Illumina Read Lengths:
    -n, --read-length NUM
          The length of the reads to
          simulate. All resulting
          reads will have the same
          length. Default: 36.

  Illumina Error Model:
    -pi, --prob-insert NUM
          Probability of an insertion.
          Default: 0.001.
    -pd, --prob-delete NUM
          Probability of a deletion.
          Default: 0.001.
    -pmmf, --prob-mismatch-file STR
          Mismatch probability path.
          If set, probability
          distribution is loaded from
          argument.
    -pmms, --prob-mismatch-scale NUM
          Scale to apply for
          probability mismatch.
          Default: 1.0.
    -pmm, --prob-mismatch NUM
          Average mismatch
          probability. Default: 0.004.
    -pmmb, --prob-mismatch-begin NUM
          Probability of a mismatch at
          the first base. Default:
          0.002.
    -pmme, --prob-mismatch-end NUM
          Probability of a mismatch at
          the last base. Default:
          0.012.
    -pr, --position-raise NUM
          Relative position of raise
          point. Default: 0.66.
    -nN, --no-N
          If set then no Ns will be
          introduced in the reads.
    -qmb, --quality-mean-begin NUM
          Quality mean at first base.
          Default: 40.
    -qme, --quality-mean-end NUM
          Quality mean at last base.
          Default: 39.5.
    -qsdb, --quality-std-dev-begin NUM
          Quality standard deviation
          at first base. Default:
          0.05.
    -qsde, --quality-std-dev-end NUM
          Quality standard deviation
          at last base. Default: 10.
    -mmqmb, --mismatch-quality-mean-begin NUM
          Mismatch quality mean at
          first base. Default: 39.5.
    -mmqme, --mismatch-quality-mean-end NUM
          Mismatch quality mean at
          last base. Default: 30.
    -mmqsdb, --mismatch-quality-std-dev-begin NUM
          Mismatch quality standard
          deviation at first base.
          Default: 3.
    -mmqsde, --mismatch-quality-std-dev-end NUM
          Mismatch quality standard
          deviation at last base.
          Default: 15.

READ NAMING SCHEME
    Note that the suffixes will not
    appear in the SAM file. In the SAM
    format, the FLAG field is used to
    mark the leftmost/rightmost
    fragment/read.

    0     No suffix, left-end and
          right-end read will be named
          'example.fastq.0000', for
          example.
    1     Add zero-based slash-suffix,
          i.e. names will end on '/0'
          and '/1'
    2     Add one-based slash-suffix,
          i.e. names will end on '/1'
          and '/2'

VERSION
    mason version: 0.1.1
    Last update July 27, 2012
