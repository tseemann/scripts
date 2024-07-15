#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use List::Util qw(max min);
use Data::Dumper;
use Math::Random qw(random_normal);

my(@Options, $verbose, $coverage, $qstart, $qend, 
             $length, $imean, $isdev, 
	     $revcom, $fastq);
setOptions();

#----------------------------------------------------------------------
# set up some static tables to save re-computation
# THIS IS A MESS!

my $pstart = 1-phred2prob($qstart);
my $pend = 1-phred2prob($qend);
my $range = abs($pend - $pstart);
my @perror = map { error_func($_/$length) } (0 .. $length-1);
@perror = map { $pend+$range*($_-$perror[-1])/($perror[0]-$perror[-1]) } @perror;

my $qstring = join('', map { chr(64-10*log(1-$perror[$_])/log(10)) } (0 .. $length-1) );
#print Dumper($pstart, $pend, \@perror, $qstring); exit;

#----------------------------------------------------------------------

my $n=0;

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
while (my $seq = $in->next_seq) {
  print STDERR "Processing ",$seq->id," ",$seq->length, " bp.\n";
  illuminate($seq);
} 
print STDERR "\n";

#----------------------------------------------------------------------

sub illuminate {
  my($seq) = @_;
  my $L = $seq->length;
  my $nread = int $L * $coverage / $length / 2;
  my @dna = ($seq->seq, $seq->revcom->seq);
  print STDERR "Generating $nread read pairs.\n";
  while ($n++ < $nread) {
    my $dir = rand() > 0.5 ? 1 : 0;
    my $start = int(rand($L - $imean - 2*$length));
    my $left = substr($dna[$dir], $start, $length);
    my $offset = random_normal(1, $imean, $isdev);
    my $begin = $start+$offset-$length;
#    die "begin=$begin L=$length Ldna=$L" if $begin > $L;
    if ($begin+$length > $L) {
      $n--;
      next;
    }
    my $right = substr($dna[$dir], $begin, $length);
    $left = mutate($left);
    $right = mutate($right);
    $right = revcom($right);
    
    next unless $left and $right;
    
    if ($revcom) { # mate pairs
      $left = revcom($left);
      $right = revcom($right); 
    }
    
    if ($fastq) {
      print "\@$imean.$n/1\n$left\n+\n$qstring\n";
      print "\@$imean.n/2\n$right\n+\n$qstring\n";
    }
    else {
      print ">$imean.$n/1\n$left\n";
      print ">$imean.$n/2\n$right\n";
    }
    
    print STDERR "\rGenerated: $n" if ($n % 9657)==0;
  }
}

#----------------------------------------------------------------------

sub mutate {
  my $s = shift;
  my $L = length $s;
  for my $i (0 .. $L-1) {
    if (rand > $perror[$i]) {
      substr($s, $i, 1) = (qw(A G T C))[int rand(4)];
    }
  }
  return $s;
}

#----------------------------------------------------------------------

sub revcom {
  my $s = shift;
  $s = reverse $s;
  $s =~ tr/ACGTacgt/TGCAtgca/;
  return $s;
} 

#----------------------------------------------------------------------
# error curve, input x is [0,1-epsilon] and output is probability [0,1]
# this approximates the illumina error curve 

sub error_func {
  my $x = shift;
  return (1 - $x) ** 0.1;
}


#----------------------------------------------------------------------

sub phred2prob {
  my $Q = shift;
  return 10**(-$Q/10);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
#    {OPT=>"n|num=i",  VAR=>\$num, DEFAULT=>36, DESC=>"Number of reads"},
    {OPT=>"c|coverage=i",  VAR=>\$coverage, DEFAULT=>50, DESC=>"Coverage"},
    {OPT=>"qstart=i",  VAR=>\$qstart, DEFAULT=>40, DESC=>"Phred quality at start of read"},
    {OPT=>"qend=i",  VAR=>\$qend, DEFAULT=>20, DESC=>"Phred quality at end of read"},
    {OPT=>"l|length=i",  VAR=>\$length, DEFAULT=>36, DESC=>"Read length"},
    {OPT=>"i|insert=i",  VAR=>\$imean, DEFAULT=>250, DESC=>"Insert size mean"},
    {OPT=>"s|sdev=f",  VAR=>\$isdev, DEFAULT=>20, DESC=>"Insert size std. dev."},
    {OPT=>"r|revcom!",  VAR=>\$revcom, DEFAULT=>0, DESC=>"Revcom reads, so opp-out for mate pair"},
    {OPT=>"q|fastq!",  VAR=>\$fastq, DEFAULT=>0, DESC=>"Output FASTQ, not FASTA"},
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
  print "Usage: $0 [options] reference.fasta > fake_illumina_reads.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
