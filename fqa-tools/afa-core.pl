#!/usr/bin/env perl
use warnings;
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Set::IntSpan;
use Data::Dumper;

my(@Options, $verbose, $informat, $outformat, $remove, $keep, $all,
             $out_fn, $core_fn, $noncore_fn, $prefix);
setOptions();

if ($prefix) {
  msg("Auto-setting all output filenames to $prefix.*");
  $out_fn = "$prefix.core.aln";
  $core_fn = "$prefix.core-sites.txt";
  $noncore_fn = "$prefix.noncore-sites.txt";
}

my $in_fn = @ARGV ? join(' ', @ARGV) : 'STDIN';
msg("Loading: $in_fn");
my $in = Bio::AlignIO->new(-fh=>\*ARGV, -format=>$informat);
my $aln = $in->next_aln;
#msg("Loaded alignment:", $aln->id);
$aln->is_flush or err("Alignment is not flush! Exiting.");
aln_stats($aln);

my $L = $aln->length;
my $N = $aln->num_sequences;

msg("Extracting $N sequences.");
my @seq;
for my $seq ($aln->each_seq) {
  push @seq, uc($seq->seq);
}

msg("Examining $L x $N-columns for matches to regexp /$remove/");
my $regexp = qr/$remove/i; # compile the regexp
my @noncore;
my $noncore=0;
COLUMN:
for my $i (0 .. $L-1) { 
  my $col = '';
  ROW:
  for my $j (0 .. $N-1) { 
    my $c = substr $seq[$j], $i, 1;
    if ($c =~ $regexp) {
      $noncore[$i] = 1;
      $noncore++;
      next COLUMN;
    }
  }
  print STDERR "\r", progress_bar($i, $L) if $i % 1000 == 0;
}
print STDERR "\n";
msg("Determined $noncore of $L columns are non-core");
msg("Core is", sprintf("%0.1f%%", ($L-$noncore)*100/$L) );

my @reject;
msg("Writing: $core_fn");
open CORE, '>', $core_fn;
msg("Writing: $noncore_fn");
open NONCORE, '>', $noncore_fn;
for my $i (0 .. $L-1) {
  if (defined $noncore[$i]) {
    print NONCORE $i+1,"\n";  # clonalframe use 1-base
    push @reject, $i;         # bioperl use 0-base
  }
  else {
    print CORE $i+1, "\n";
  }
  print STDERR "\r", progress_bar($i, $L) if $i % 1000 == 0;
}
print STDERR "\n";

my $core_aln;
if ($noncore > 0) {
  msg("Removing $noncore columns from alignment. Please be patient...");
  my $set = Set::IntSpan->new(@reject);
  #print Dumper($set->spans);
  $core_aln = $aln->remove_columns( $set->spans );
}
else {
  msg("No non-core sites found!");
  $core_aln = $aln;
}
msg("Core alignment has", $core_aln->length, "columns");

msg("Restoring original sequence names.");
$core_aln->set_displayname_flat();

msg("Writing core alignment: $out_fn");
my $out = Bio::AlignIO->new(-file=>">$out_fn", -format=>$outformat);
$out->write_aln($core_aln);

msg("Done.");

#.......................................................................
sub progress_bar {
  my($count, $total, $prefix, $width) = @_;
  $prefix ||= "Progress: ";
  $width = $width || $ENV{COLUMNS} || 70;
  $width = $width - length($prefix) - 2;  # 2 = [ ..... ]
  my $done = int($count * $width / $total + 0.5);
  my $bar = '='x($done) . '.'x($width-$done);
  return "$prefix\[$bar\]";
}

#.......................................................................
sub aln_stats {
  my($aln) = @_;
  msg("Sequences:", $aln->num_sequences);
  msg("Length:", $aln->length);
#  printf STDERR "Non-gaps: %d\n", $aln->num_residues;
#  printf STDERR "Alphabet: %s\n", join(' ', $aln->symbol_chars);
}

#.......................................................................
sub msg {
  my $t = localtime;
  print STDERR "@_\n";
}
            
#.......................................................................
sub err {
  msg(@_);
  exit(1);
}
                
#.......................................................................
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"if=s",  VAR=>\$informat, DEFAULT=>'fasta', DESC=>"Input format"},
    {OPT=>"of=s",  VAR=>\$outformat, DEFAULT=>'fasta', DESC=>"Output format"},
    {OPT=>"remove=s",  VAR=>\$remove, DEFAULT=>'[^AGTC]', DESC=>"Removes columns matching this regexp"},
#    {OPT=>"keep=s",  VAR=>\$keep, DEFAULT=>'AGTC', DESC=>"Keep columns with these characters"},
#    {OPT=>"all!",  VAR=>\$all, DEFAULT=>0, DESC=>"Require column to contain *only* those chars"},
    {OPT=>"prefix=s",  VAR=>\$prefix, DEFAULT=>'', DESC=>"Output all the files with this name"},
    {OPT=>"out=s",  VAR=>\$out_fn, DEFAULT=>'/dev/stdout', DESC=>"Output alignment file"},
    {OPT=>"core=s",  VAR=>\$core_fn, DEFAULT=>'/dev/null', DESC=>"Output file core sites (for ClonalFrameML)"},
    {OPT=>"noncore=s",  VAR=>\$noncore_fn, DEFAULT=>'/dev/null', DESC=>"Output file non-core sites (ClonalFrameML)"},
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
