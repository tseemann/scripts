#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use Bio::Seq::Quality;

my(@Options, $verbose, $library, $prefix, $append, $nolib);
setOptions();

# HOW TO:
#
# phred -id chromat_dir/ -pd phd_dir/ -log
# rm -f 454Sanger.fna 454Sanger.fna.qual
# find phd_dir -name \*phd.1 -print | sort | xargs phred-to-454sanger.pl --append

# NEWBLER FASTA SPECS:
#
#1. template - the Paired End template string for this read (Paired End reads are matched by
#   having the same template string)
#2. dir - values "F", "R", "fwd" or "rev" giving the direction of the Paired End read
#3. library - the name of the library that generated this Paired End read (all Paired End reads
#   are grouped by library name for the determination of expected pair distance)
#4. trim - the trimmed region of the sequence, given as "#-#"
#5. scf - the path or "command string" to use to access the SCF file for the read
#6. phd - the path or "command string" to use to access the PHD file for the read

$prefix or die "please supply --prefix PREFIX";
$library = $prefix.scalar(@ARGV) unless $library;
print STDERR "Setting library=$library\n";

my $mode = $append ? '>>' : '>';
print STDERR "Writing ($mode): $prefix.fna $prefix.fna.qual\n";
my $so = Bio::SeqIO->new(-file=>"$mode$prefix.fna", -format=>'fasta');
my $qo = Bio::SeqIO->new(-file=>"$mode$prefix.fna.qual", -format=>'qual');
my $count=0;
my $TOTAL = scalar(@ARGV);

for my $file (@ARGV) {
  my $phdio = Bio::SeqIO->new(-file=>$file, -format=>'phd');
  my $phd = $phdio->next_seq;
#  print Dumper($phd);
  my $trim = $phd->attribute('TRIM');
#  print Dumper($trim);
  my($begin,$end) = split ' ', $trim;
  next if $begin < 0 or $end < 0 or $phd->length < 1;
  $count++;
  my $id = $phd->display_id;
  #
  # this will have to be customised for every set of chromats you get!!
  $id =~ m/^(.*?)(\.s.\d)?\.([bg]\d*)$/;
  #
  my $template = $1;
  my $dir = ($3 =~ m/[bwf]/ ? 'fwd' : 'rev');
  my $lib = $nolib ? '' : "library=$library ";
  $phd->desc("${lib}template=$template dir=$dir trim=$begin-$end");
  print STDERR "\r[$count/$TOTAL] $id ", $phd->desc;
  $so->write_seq($phd); # dna
  $qo->write_seq($phd); # quality
}
print STDERR "\nWrote $TOTAL sequences.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"p|prefix=s", VAR=>\$prefix, DEFAULT=>'454Sanger', DESC=>"Output file prefix"},
    {OPT=>"a|append=s", VAR=>\$append, DEFAULT=>0, DESC=>"Append to PREFIX files"},
    {OPT=>"l|library=s", VAR=>\$library, DEFAULT=>'', DESC=>"Library name"},
    {OPT=>"n|nolib!", VAR=>\$nolib, DEFAULT=>0, DESC=>"Sequences not from library eg. PCR products"},
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
  print "Usage: $0 [options] trace1.phd [trace2.phd ...]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
