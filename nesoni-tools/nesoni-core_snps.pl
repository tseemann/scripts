#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use List::MoreUtils qw(pairwise);
use List::Util qw(sum);

no warnings qw(once); #Name "main::a" used only once: possible typo

my(@Options, $verbose);
setOptions();

my $refname = "$ARGV[0]/reference.fa";
print STDERR "Loading reference $refname\n";
my $refio = Bio::SeqIO->new(-file=>$refname, -format=>'fasta');
my $ref = $refio->next_seq->seq;
my $len = length($ref);
print STDERR "Loaded $len bp.\n";

# initialise with 1s - we will multiply by coverage at each point
# any zero coverage will zero it forever
my @cov = map { 1 } (1 .. $len);

for my $ndir (@ARGV) {
  next unless -d $ndir;
  my($fn) = glob("$ndir/*-ambiguous-depth.userplot");
  print STDERR "Loading $fn - ";
  if (not -r $fn) {
    print STDERR " could not read, skipping.\n";
  }
  else {
    open my $fh, '<', $fn;
    my @depth = <$fh>;
    close $fh;
    printf STDERR "%d lines\n", scalar @depth;
    if (@depth != $len) {
      print STDERR "ERROR: $fn does not have length $len bp - skipping\n";
      next;
    }
    # the trick!
    chomp @depth;
    @cov = pairwise { $a * ($b > 0) } @cov, @depth;
    print STDERR "Processed $ndir, core is now ",sum(@cov),"/$len\n";
  }
}

printf ">CoreGenome $len bp\n";
my $ncore=0;

for my $i (0 .. $len-1) {
  print STDERR "\rWriting: $i/$len bp" if not $verbose and $i%971==0;
  print $cov[$i] ? do { $ncore++; substr($ref,$i,1) } : 'N';
  print "\n" if (($i+1)%60==0);
}
print "\n";
printf STDERR "\nCore genome is $ncore/$len bases (%.1f%%)\n", $ncore*100/$len;


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
  print "Usage: $0 [options] nesoniDir1 [nesoniDir2 ...] > core_genome.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
