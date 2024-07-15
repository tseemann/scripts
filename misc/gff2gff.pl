#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::Tools::GFF;

my(@Options, $verbose, $infile, $outfile, $informat, $outformat);
setOptions();

my $in = Bio::Tools::GFF->new(-file=>$infile, -gff_version=>$informat);
my $out = Bio::Tools::GFF->new(-file=>">$outfile", -gff_version=>$outformat);

$in->ignore_sequence(0); 
$in->features_attached_to_seqs(1);

while (my $f = $in->next_feature) {
  print STDERR "\rProcessing feature: $.";
  $out->write_feature($f);
}
print STDERR "\n";

for my $seq ($in->get_seqs) {
  print STDERR "\rProcessing sequence: $.";
  $out->write_seq($seq);
}
print STDERR "\nDone.\n";


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"i|infile=s", VAR=>\$infile, DEFAULT=>'/dev/stdin', DESC=>"Input filename"},
    {OPT=>"o|outfile=s", VAR=>\$outfile, DEFAULT=>'/dev/stdout', DESC=>"Output filename"},
    {OPT=>"if|informat=i", VAR=>\$informat, DEFAULT=>2, DESC=>"Input GFF format: 1,2,3 [use 2 for GTF]"},
    {OPT=>"of|outformat=i", VAR=>\$outformat, DEFAULT=>3, DESC=>"Output GFF format: 1 2 3 [use 2 for GTF]"},
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
  print "Usage: $0 [options] -i file.gtf -if 2 -o file.gff3 -of 3\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

