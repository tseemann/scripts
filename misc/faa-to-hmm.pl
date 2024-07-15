#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use Bio::AlignIO;
use File::Temp qw(tempdir);
use File::Copy;

my(@Options, $verbose, $aln_id, $aln_acc, $aln_desc, $outname);
setOptions();

my $dir = tempdir(CLEANUP=>1);

my $nseq=0;
my $fasta = "$dir/seq.faa";
print STDERR "Preparing sequences: $fasta\n";
my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-file=>">$fasta", -format=>'Fasta');
while (my $seq = $in->next_seq) {
  $out->write_seq($seq);
  $nseq++;
} 
print STDERR "Found sequences: $nseq\n";

my $cmd = "muscle -quiet -in $fasta";
print STDERR "Aligning sequences: $cmd\n";
my $muscle = Bio::AlignIO->new(-file=>"$cmd |", -format=>'fasta');
my $aln = $muscle->next_aln;

$aln->id($aln_id) if $aln_id;
$aln->accession($aln_acc) if $aln_acc;
$aln->description($aln_desc) if $aln_desc;

my $stockholm = "$dir/seq.aln";
print STDERR "Writing stockholm alignment: $stockholm\n";
my $aout = Bio::AlignIO->new(-file=>">$stockholm", -format=>'stockholm');
$aout->write_aln($aln);

my $hmmer = "$dir/seq.hmm";
$cmd = "hmmbuild -o /dev/null $hmmer $stockholm";
print STDERR "Building HMM: $cmd\n";
system($cmd);

if ($outname) {
  print STDERR "Writing: $outname\n";
  move($hmmer, $outname);
}
else {
  system("cat $hmmer");
}
print STDERR "Cleaning up.\nDone.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"id=s",  VAR=>\$aln_id, DEFAULT=>'', DESC=>"Set HMM id"},
    {OPT=>"acc=s",  VAR=>\$aln_acc, DEFAULT=>'', DESC=>"Set HMM accession"},
    {OPT=>"desc=s",  VAR=>\$aln_desc, DEFAULT=>'', DESC=>"Set HMM description"},
    {OPT=>"o|out=s",  VAR=>\$outname, DEFAULT=>'', DESC=>"Set output file instead of STDOUT"},
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
