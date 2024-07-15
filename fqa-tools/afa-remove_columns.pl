#!/usr/bin/env perl
use warnings;
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;

my(@Options, $verbose, $informat, $outformat, $char, $all);
setOptions();

$char ||= 'N';
$char = substr $char, 0, 1;

my $in  = Bio::AlignIO->new(-fh=>\*ARGV, -format=>$informat);
my $out = Bio::AlignIO->new(-fh=>\*STDOUT, -format=>$outformat);

my $aln = $in->next_aln;
printf STDERR "Loaded alignment: %s\n", $aln->id;
aln_stats($aln);

my $constraint = $all ? 'all/only' : 'at least one';
print STDERR "Removing columns containing $constraint '$char' ...\n";
#$aln->gap_char($char);
$aln = $aln->remove_gaps($char, $all);
aln_stats($aln);

print STDERR "Restoring original sequence names...\n";
$aln->set_displayname_flat();

print STDERR "Writing new alignment...\n";
$out->write_aln($aln);

print STDERR "Done.\n";

#.......................................................................

sub aln_stats {
  my($aln) = @_;
  printf STDERR "Sequences: %d\n", $aln->num_sequences;
  printf STDERR "Length: %d\n", $aln->length;
#  printf STDERR "Non-gaps: %d\n", $aln->num_residues;
  printf STDERR "Alphabet: %s\n", join(' ', $aln->symbol_chars);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"if=s",  VAR=>\$informat, DEFAULT=>'fasta', DESC=>"Input format"},
    {OPT=>"of=s",  VAR=>\$outformat, DEFAULT=>'fasta', DESC=>"Output format"},
    {OPT=>"char=s",  VAR=>\$char, DEFAULT=>'N', DESC=>"Removes columns with this character"},
    {OPT=>"all!",  VAR=>\$all, DEFAULT=>0, DESC=>"Require column to contain *only* --char symbols'"},
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
