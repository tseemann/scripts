#!/usr/bin/env perl
use strict;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq assert int_qual write_fasta);

my(@Options, $verbose, $sfile, $qfile, $offset);
setOptions();

print STDERR "Opening: @ARGV\n";

assert($sfile, "Please supply --fasta");
open my $sfh, '>', $sfile;
print STDERR "Writing: $sfile\n";

$qfile = $sfile.".qual" if $qfile eq 'auto';
open my $qfh, '>', $qfile;
print STDERR "Writing: $qfile\n";

while ( not eof() ) {
  my $r = read_fastq(\*ARGV);
  write_fasta($sfh, $r);
  my @q = int_qual($offset, $r);
  @q = map { $_ < 0 ? 0 : $_ } @q;  # make -ve ones zero
  print $qfh ">",$r->[0],"\n@q\n";
}
print STDERR "Done.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"f|fasta=s", VAR=>\$sfile, DEFAULT=>'', DESC=>"Output FASTA file"},
    {OPT=>"q|qual=s", VAR=>\$qfile, DEFAULT=>'auto', DESC=>"Output quality file"},
    {OPT=>"o|offset=i", VAR=>\$offset, DEFAULT=>32, DESC=>"FASTQ quality offset: sanger=32 illumina=64"},
  );

  (@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 input.fq -f output.fasta -q output.qual\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
