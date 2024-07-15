#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use File::Path;
use Time::Piece;
use Cwd;

# Perl Wrapper to ease pain of using Minimus2 AMOS/nucmer-based assembler
#
# http://sourceforge.net/apps/mediawiki/amos/index.php?title=Minimus2
#
#cat 454Scaffolds.fna contigs.fa > mini_in_orig_scaffs.fa
#toAmos -s mini_in_orig_scaffs.fa -o mini_in_orig_scaffs.afg
#minimus2 mini_in_orig_scaffs -D REFCOUNT=25
#
#minimus2 prefix  \
#   -D REFCOUNT=n  \  # Number of sequences is the first set
#   -D OVERLAP=n   \  # Minimum overlap (Default 40bp)
#   -D CONSERR=f   \  # Maximum consensus error (0..1) (Def 0.06)
#   -D MINID=n     \  # Minimum overlap %id for align. (Def 94)
#   -D MAXTRIM=n      # Maximum sequence trimming length (Def 20bp)


my(@Options, $dir, $clean, $verbose, $overlap, $conserr, $minid, $maxtrim);
setOptions();

my $STEM = 'out';  # to be compatible with default nucmer stem

if (not $dir) {
  my $t = localtime;
  $dir = "minimus_".$t->datetime;
  msg("Automatic dir name: $dir");
}

msg("Creating folder: $dir");
mkpath($dir, 1);

my $FASTA = "$dir/$STEM.fa";
msg("Will write to: $FASTA");

my $allseq = Bio::SeqIO->new(-file=>">$FASTA", -format=>'Fasta');
my @numseq;
for my $index (0 .. 1) {
  msg("Reading: $ARGV[$index]");
  my $fa = Bio::SeqIO->new(-file=>$ARGV[$index], -format=>'Fasta');
  while (my $seq = $fa->next_seq) {
    $allseq->write_seq($seq);
    $numseq[$index]++;
  }
  msg("Sequences: $numseq[$index]");
} 

my $AFG = "$dir/$STEM.afg";
msg("Running toAmos to create $AFG");
system("toAmos -s \Q$FASTA\E -o \Q$AFG\E");

my $oldcwd = getcwd;
msg("Changing directory: $dir");
chdir($dir);
msg("Running minimus2");
system(
  "minimus2 $STEM".
  " -D REFCOUNT=$numseq[0]".
  " -D OVERLAP=$overlap".
  " -D CONSERR=$conserr".
  " -D MINID=$minid".
  " -D MAXTRIM=$maxtrim"
);
msg("Output files:");
system("ls -s $STEM.*");

msg("Concerning lines in: $STEM.runAmos.log");
system("egrep -i 'forced|error' $STEM.runAmos.log");

if ($clean) {
  msg("Cleaning up.");
#  rmtree("$STEM.bnk");
  unlink "$STEM.fa", "$STEM.OVL";
}

msg("Returning to dir: $oldcwd");
chdir $oldcwd;

msg("Generating result summary:");
my @look = ($ARGV[0], $ARGV[1], "$dir/$STEM.fasta", "$dir/$STEM.singletons.seq");
#print "[ALL]\n";
for my $file (@look) {
  system("fa -e \Q$file\E");
}
#print "[LARGE >= 500bp]\n";
#for my $file (@look) {
#  system("fa -e --minsize 500 \Q$file\E");
#}

msg("Please examine $dir/$STEM.runAmos.log for alignment exception information.");
msg("Done.");

sub msg {
  print STDERR @_,"\n" if $verbose;
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>1, DESC=>"Verbose output"},
    {OPT=>"d|dir=s",  VAR=>\$dir, DEFAULT=>'', DESC=>"Directory to put output (default=minimus_YYYY-MM-DDThh:mm:ss)"},
    {OPT=>"c|clean!",  VAR=>\$clean, DEFAULT=>0, DESC=>"Clean excess output from --dir"},
    {OPT=>"o|overlap=i",  VAR=>\$overlap, DEFAULT=>40, DESC=>"Minimum overlap"},
    {OPT=>"e|conserr=f",  VAR=>\$conserr, DEFAULT=>0.06, DESC=>"Maximum consensus error"},
    {OPT=>"m|minid=i",  VAR=>\$minid, DEFAULT=>94, DESC=>"Minimum overlap \%id for align."},
    {OPT=>"t|maxtrim=i",  VAR=>\$maxtrim, DEFAULT=>20, DESC=>"Maximum sequence trimming length"},
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
  print "Usage: $0 [options] [-d folder] seqs_1.fasta seqs_2.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
