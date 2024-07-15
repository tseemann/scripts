#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use File::Temp;

my(@Options, $verbose, $dir);
setOptions();

die "THIS PROGRAM IS INCOMPLETE AND BUGGY, SORRY!";

-d $dir or die "can't read chromat_dir named '$dir'";
opendir DIR, $dir;
my @abi = grep { !/^\./ && -r "$dir/$_" } readdir(DIR);
closedir DIR;

print STDERR "Found ",0+@abi," trace files.\n";

my $list = File::Temp->new(SUFFIX => '.list');
my $seq = File::Temp->new(SUFFIX => '.fasta');
my $qual = File::Temp->new(SUFFIX => '.qual');

select((select($list), $| = 1)[0]);
select((select($seq), $| = 1)[0]);
select((select($qual), $| = 1)[0]);

print STDERR "Temps: $list $seq $qual\n";

for my $f (@abi) {
  seek $list, 0, 'SEEK_SET';
  seek $seq, 0, 'SEEK_SET';
  seek $qual, 0, 'SEEK_SET';
  print $list "$dir/$f\n";
  system("phred -if $list -s $seq -q $qual");
#  exit;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"d|chromatdir=s", VAR=>\$dir, DEFAULT=>'', DESC=>"Examine files in chromat_dir"},
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
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
