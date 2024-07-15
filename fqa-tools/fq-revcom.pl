#!/usr/bin/env perl
use strict;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fastq assert reverse_complement);
use Data::Dumper;

my(@Options, $verbose, $suffix, $prefix);
setOptions();

my $THROB = 7931;

my $input = @ARGV ? "$ARGV[0] ..." : "(stdin)";
print STDERR "Reading: $input\n";

my $nread=0;
my $seq;

while ( not eof () ) {
  $seq = read_fastq(\*ARGV);
  reverse_complement($seq);
  $seq->[0] = $prefix.$seq->[0].$suffix;
  write_fastq(\*STDOUT, $seq);
  print STDERR "\rProcessing: $nread" if $nread % $THROB == 0;
  $nread++;
}

print STDERR "\rProcessed: $nread     \n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbosity"},
    {OPT=>"p|prefix=s",   VAR=>\$prefix, DEFAULT=>'', DESC=>"Prepend this prefix to the sequence ID"},
    {OPT=>"s|suffix=s",   VAR=>\$suffix, DEFAULT=>'', DESC=>"Append this prefix to the sequence ID"},
  );

  #(@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 reads.fq > revcom_reads.fq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
