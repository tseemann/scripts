#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Data::Dumper;
use List::Util qw(sum min max);

my(@Options, $debug, $ends, $flank, $minext);
setOptions();

sub read_seq {
  my($fh) = @_;
  my $s = '';
  while (<$fh>) {
    chomp;
    last if $_ eq '';
    $s .= $_;
  }
  return $s;
}
    
my $cid;
my %contig;

while (<ARGV>) {
  if (m/^CO (\S+)/) {
    $cid = $1;
    print STDERR "Reading contig: $cid\n";
    $contig{$cid}{SEQ} = read_seq(\*ARGV);
    $contig{$cid}{LEN} = length $contig{$cid}{SEQ};
  }
  elsif (m/^AF (\S+) ([UC]) (-?\d+)/) {
    $contig{$cid}{READS}{$1}{POS} = $3;
  }
  elsif (m/^RD (\S+) (\d+)/) {
    $contig{$cid}{READS}{$1}{SEQ} = read_seq(\*ARGV);
    $contig{$cid}{READS}{$1}{LEN} = length $contig{$cid}{READS}{$1}{SEQ};
  }
}
#print STDERR Dumper(\%contig) if $debug;
print STDERR "Contigs: ", scalar(keys %contig), "\n";
print STDERR "Reads: ", sum(map { scalar(keys %{$contig{$_}{READS}}) } keys %contig), "\n";

my $fasta = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $bp=0;
for my $cid (sort keys %contig) {
  my $ctg = $contig{$cid};
  my $minx = 0;
  my $maxx = 0;
  my $minid = '';
  my $maxid = '';
  for my $rid ( keys %{$ctg->{READS}} ) {
    my $read = $ctg->{READS}{$rid};
    my $pos = $read->{POS};
    my $extra = 0;
    if ($pos < 1) {
      $extra = $pos;
    }
    elsif ($pos+$read->{LEN} > $ctg->{LEN}) {
      $extra = $pos + $read->{LEN} - $ctg->{LEN};
    }
    printf STDERR "%s %+4d\n", $cid, $extra if $extra and $debug;
    if ($extra < $minx) {
      $minx = $extra;
      $minid = $rid;
    }
    elsif ($extra > $maxx) {
      $maxx = $extra;
      $maxid = $rid;
    }
#    $minx = min($extra, $minx);
#    $maxx = max($extra, $maxx);    
  }
  print STDERR "$cid\t$minid\t$minx\t$maxid\t$maxx\n" if $debug;
  $bp += abs($minx) + abs($maxx);
  
  my $left  = $minid ? substr( $ctg->{READS}{$minid}{SEQ}, 0, abs($minx) ) : '';
  my $right = $maxid ? substr( $ctg->{READS}{$maxid}{SEQ}, $ctg->{READS}{$maxid}{LEN}-$maxx, $maxx ) : '';
  print STDERR ">LEFT\n$left\n" if $debug;
  print STDERR ">RIGHT\n$right\n" if $debug;
  my $middle = $ctg->{SEQ};
  $middle =~ s/\*//g;
  if ($ends) {
    if ($left and length($left) >= $minext) {
      my $dna = lc($left);
      $dna .= uc substr($middle, 0, $flank);
      $fasta->write_seq( Bio::Seq->new( -id=>"$cid.ext.L", -seq=>$dna ) );
    }
    if ($right and length($right) >= $minext) {
      my $dna = lc($right);
      $dna = (uc substr($middle, -$flank, $flank)) . $dna;
      $fasta->write_seq( Bio::Seq->new( -id=>"$cid.ext.R", -seq=>$dna ) );
    }
  }
  else {
    $fasta->write_seq(
      Bio::Seq->new( -id=>"$cid.ext", -seq=>lc($left).uc($middle).lc($right) )
    );
  }
}

print STDERR "Could 'reclaim' $bp bp.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"flank=i",  VAR=>\$flank, DEFAULT=>100, DESC=>"Add this much flanking contig sequence in --ends mode"},
    {OPT=>"minext=i",  VAR=>\$minext, DEFAULT=>100, DESC=>"Minimum length of extension to consider"},
    {OPT=>"ends!",  VAR=>\$ends, DEFAULT=>0, DESC=>"Output just ends, not whole new contigs"},
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
  print "Usage: $0 454Contigs.ace > 454ContigsReclaimed.fasta\n";
  print "       $0 --ends 454Contigs.ace > 454ContigsMissingEnds.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
