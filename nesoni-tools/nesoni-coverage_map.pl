#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;

my(@Options, $verbose, $purity, $depth, $folder, $trim);
setOptions();

for my $file (@ARGV) {
  problem("can't read file '$file'") unless -r $file;
}

my $ref = shift @ARGV;
msg("Reference file: $ref");

my @reads = @ARGV;
msg("Read file(s): @reads");

if (not $folder) {
  $folder = "nesoni-$ref";
}

msg("Output folder: $folder/");
problem("output folder '$folder' already exists") if -d $folder;

my $cmd = "nice nesoni shrimp $folder $ref".
          " reads @reads shrimp-options -R";
runcmd($cmd);

$cmd = 
  "nesoni consensus $folder --depth $depth --purity $purity".
  " --ambiguity-codes 1 --fidelity 1.00 --monogamous yes".
  " --max-pair-sep 500".
  "";
runcmd($cmd);
my($snps) = qx(cat "$folder/report.txt" | wc -l);
chomp $snps;
$snps--; # header line
msg("Found $snps SNPs in $folder/report.txt");

my($firstline) = qx(head -1 $ref);
if ($firstline =~ m/^LOCUS/) {
  msg("$ref looks like a Genbank file, running consequences");
  $cmd = "nesoni consequences --tabular --transl_table 11 $ref $folder";
  runcmd($cmd);
}
else {
  msg("$ref does not look like a Genbank file");
}

#----------------------------------------------------------------------

sub runcmd {
  my(@cmd) = @_;
  print STDERR "Running: @cmd\n";
  system(@cmd)==0 or problem("system(@cmd) failed: $?");
}

#----------------------------------------------------------------------

sub problem {
  print STDERR "ERROR:\n@_\n";
  exit -1;
}

#----------------------------------------------------------------------

sub msg {
  print STDERR "@_\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"d|depth=i",  VAR=>\$depth, DEFAULT=>2, DESC=>"Nesoni consensus --depth"},
    {OPT=>"p|purity=f",  VAR=>\$purity, DEFAULT=>0.5, DESC=>"Nesoni consensus --purity"},
    {OPT=>"t|trim=i",  VAR=>\$trim, DEFAULT=>2, DESC=>"Nesoni consensus --trim"},
    {OPT=>"f|folder=s",  VAR=>\$folder, DEFAULT=>'', DESC=>"Destination folder (default = <REFERENCE_FILENAME>.nesoni)"},
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
  print "Usage: $0 [options] <ref.fna | ref.gbk>  <reads.fa | reads.fq>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
