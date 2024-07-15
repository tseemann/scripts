#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;

# SAMSHRIMP
#
#    --cs yes/no            - Colorspace reads (i.e. use gmapper-cs).
#                             Default: no
#    --sam-unaligned yes/no - Pass --sam-unaligned to gmapper
#                             Default: yes   
#    --half-paired yes/no   - Pass --half-paired to gmapper
#                             Default: yes

my(@Options, $verbose, $folder, $solid, $ungapped, $single, $cutoff, 
             $exe, $unmapped, $pcid, $whole, $strict);
setOptions();

if ($pcid < 1 or $pcid > 100) {
  problem("--pcid must be between 1 and 100");
}

for my $file (@ARGV) {
  problem("can't read file '$file'") unless -r $file;
}

my $ref = shift @ARGV or problem("No reference supplied");
msg("Reference file: $ref");

my @reads = @ARGV or problem("No reads supplied"); 
msg("Read file(s): @reads");

if (not $folder) {
  my $name = $ref;
  $name =~ s{^(.*/)}{};
  $folder = "nesoni-$ref";
}

if (!$single and $reads[0] =~ m/single/i) {
  $single=1;
  msg("*** Turning --single ON, as $reads[0] has 'single' in name!");
}

if ($single and $reads[0] =~ m/pair/i) {
  $single=0;
  msg("*** Turning --single OFF, as $reads[0] has 'pair' in name!");
}

msg("Output folder: $folder/");
problem("output folder '$folder' already exists") if -d $folder;

if ($solid) {
  msg("Enabling SOLiD support.");
  $solid = "--solid";
}
else {
  $solid = '';
}

if ($unmapped) {
  $unmapped = "--sam-unaligned yes";
}
else {
  $unmapped = "--sam-unaligned no";
}

my $inmode = $single ? 'reads' : 'interleaved';
my $cmd = "nice $exe samshrimp $folder $unmapped $ref $inmode @reads";
my @opt;
push @opt, "--ungapped" if $ungapped;
push @opt, "-p opp-in -I -50,500" unless $single;
push @opt, "-h $pcid\%";
push @opt, "--half-paired no" if $strict;
$cmd .= " shrimp-options @opt" if @opt;
runcmd($cmd);

$whole = $whole ? "--whole-read-only yes" : "";

@opt = ("--ambiguity-codes no");
push @opt, "--majority 0.90" if $strict;
push @opt, "--whole-read-only yes" if $strict or $whole;
$cmd = "$exe samconsensus @opt $folder";
runcmd($cmd);
my($snps) = qx(cat "$folder/report.txt" | wc -l);
chomp $snps;
$snps--; # header line
msg("Found $snps SNPs in $folder/report.txt");

my($firstline) = qx(head -1 $ref);
if ($firstline =~ m/^LOCUS/) {
  msg("$ref looks like a Genbank file, running consequences");
  $cmd = "nesoni consequences --tabular $ref $folder > $folder/consequences.csv"; # pypy will fail
  runcmd($cmd);
}
else {
  msg("$ref does not look like a Genbank file, so we didn't run consequences on it.");
}

#----------------------------------------------------------------------

sub runcmd {
  my(@cmd) = @_;
  print STDERR "Running: @cmd\n";
  system(@cmd)==0 or problem("system(@cmd) failed: $?");
}

#----------------------------------------------------------------------

sub problem {
  print STDERR "ERROR: @_\n";
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
    {OPT=>"u|ungapped!",  VAR=>\$ungapped, DEFAULT=>0, DESC=>"Don't look for indels, only subs (runs faster)"},
    {OPT=>"1|single!",  VAR=>\$single, DEFAULT=>0, DESC=>"Reads are SE, not PE (will auto-detect too)"},
    {OPT=>"f|folder=s",  VAR=>\$folder, DEFAULT=>'', DESC=>"Destination folder (default = nesoni-<REFERENCE_FILENAME>)"},
    {OPT=>"e|exe=s",  VAR=>\$exe, DEFAULT=>'pypy-nesoni', DESC=>"Nesoni command: pypy-nesoni OR nesoni"},
    {OPT=>"n|unmapped!",  VAR=>\$unmapped, DEFAULT=>0, DESC=>"Include unmapped reads"},
    {OPT=>"p|pcid=i",  VAR=>\$pcid, DEFAULT=>68, DESC=>"Percent of max score needed to count"},
    {OPT=>"w|wholereads!",  VAR=>\$whole, DEFAULT=>0, DESC=>"Force --whole-read-only"},
    {OPT=>"s|strict!",  VAR=>\$strict, DEFAULT=>0, DESC=>"Strict mode: --pcid 90 --wholereads --nohalfpaired --majority 0.9"},
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
