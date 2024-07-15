#!/usr/bin/env perl
use strict;
use Data::Dumper;

my(@Options, $verbose, $prefix, $ungapped, $single, $cpus, $insert);
setOptions();

my $ref = shift @ARGV or problem("No reference supplied");
msg("Reference file: $ref");

my $reads = shift @ARGV or problem("No read file supplied"); 
msg("Read file: $reads");

if (not $prefix) {
  $ref =~ m{^(.*?/)?(\S+?)\..*$};
  $prefix = $2 || 'gmapper';
}
msg("Using prefix: $prefix");

msg("Mapping reads: $prefix.sam");
$ungapped = $ungapped ? '-U' : '';
$single = $single ? '' : "-p opp-in -I 0,$insert";
my $cmd = "nice gmapper -N $cpus -E $single $ungapped $reads $ref > $prefix.sam 2> /dev/null";
runcmd($cmd);

msg("Converting SAM to BAM: $prefix.unsorted.bam");
$cmd = "samtools view -b -t $ref -o $prefix.unsorted.bam $prefix.sam";
runcmd($cmd);

msg("Sorting BAM: $prefix.bam");
# it adds the suffix .bam for you
$cmd = "samtools sort $prefix.unsorted.bam $prefix"; 
runcmd($cmd);

msg("Indexing BAM: $prefix.bam.*");
$cmd = "samtools index $prefix.bam";
runcmd($cmd);

for my $delme ("$prefix.sam", "$prefix.unsorted.bam") {
  msg("Deleting temporary file: $delme");
  unlink $delme;
}


#         -->   -->
#         GGA   GTA
#      + TGGATCCGTAA
#      - ACCTAGGCATT
#         CCT   CAT
#         <--   <--
#    (R1:GGA, R2:TAC) or (R1:TAC, R2:GGA) are "opp-in"
#    (R1:GTA, R2:TCC) or (R1:TCC, R2:GTA) are "opp-out"
#    (R1:GGA, R2:GTA) or (R1:TAC, R2:TCC) are "col-fw"
#    (R1:GTA, R2:GGA) or (R1:TCC, R2:TAC) are "col-bw"

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
    {OPT=>"p|prefix=s",  VAR=>\$prefix, DEFAULT=>'', DESC=>"Out file prefix: .bam .bam.idx ..."},
    {OPT=>"s|single!",  VAR=>\$single, DEFAULT=>0, DESC=>"Assume single-end reads"},
    {OPT=>"c|cpus=i",  VAR=>\$cpus, DEFAULT=>8, DESC=>"Use this many CPUs"},
    {OPT=>"i|insert=i",  VAR=>\$insert, DEFAULT=>500, DESC=>"Insert size for paired reads"},
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
  print "Usage: $0 [options] ref.fasta pairs.fa          # paired reads\n";
  print "       $0 [options] --single ref.fasta reads.fa # fragment library\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
