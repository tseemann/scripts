#!/usr/bin/env perl
use warnings;
use strict;

my(@Options, $verbose, $upto);
setOptions();

#Col Field Type Regexp/Range Brief description
#1 QNAME String [!-?A-~]f1,255g Query template NAME
#2 FLAG Int [0,2^16-1] bitwise FLAG
#3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
#4 POS Int [0,2^29-1] 1-based leftmost mapping POSition
#5 MAPQ Int [0,2^8-1] MAPping Quality
#6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
#7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next segment
#8 PNEXT Int [0,229-1] Position of the mate/next segment
#9 TLEN Int [-229+1,229-1] observed Template LENgth
#10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
#11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33

my %len_of;
my %pair;
my $unmapped=0;
my $mapped=0;
my $seen=0;
my $numrev=0;

while (<ARGV>) {
  if ( m/^\@SQ\tSN:(.+)\tLN:(\d+)/ ) {   # bowtie2 order...
    $len_of{$1} = $2;
    next;
  }
  next if m/^\@/;
  my($id,$flags,$contig,$pos) = split m/\t/;
#  next unless exists $len_of{$contig};
  if ($contig eq '*') {
    $unmapped++;
  }
  else {
#    print "$id hit $contig at $pos (unmapped=$unmapped)\n";
    $mapped++;
    my $rev = $flags & 0x10;
    $numrev++ if $rev;
  }
  $seen++;
  last if $upto and ++$seen >= $upto;
  print STDERR "\rProcessed: $seen" if $seen % 10_000 == 0;
}

print STDERR "\n";
print "Read $seen, mapped $mapped ($numrev rev), unmapped $unmapped\n";


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"upto=i",  VAR=>\$upto, DEFAULT=>0, DESC=>"Only process this many lines (0=all)"},
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
  print "Usage: $0 [options] <input.sam>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
