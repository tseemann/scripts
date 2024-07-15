#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Data::Dumper;
use IO::File;

my(@Options, $verbose, $desc, $bases, $full, $minsize, $each, $width, $tabbed);
setOptions();

$full=1 if $bases;

my %inputs = $each ? (map { ($_ => IO::File->new($_)) } @ARGV) : ( '(stdin)' => \*ARGV ) ;

print "Name\tno\tbp\tok\tNs\tgaps\tmin\tavg\tmax\tN50\n" if $tabbed;

for my $fname (sort keys %inputs) 
{
  my $numseq=0;
  my $avglen=0;
  my $minlen=1E9;
  my $maxlen=0;
  my @len;
  my $toosmall=0;
  my $nn=0;
  my $xx=0;

#  print STDERR "Reading: $fname\n" if $verbose;
#  my $in = Bio::SeqIO->new(-file=>$fname, -format=>'Fasta') or next;
  my $in = Bio::SeqIO->new(-fh=>$inputs{$fname}, -format=>'Fasta') or next;
  while (my $seq = $in->next_seq) {
    my $L = $seq->length;
    if ($L < $minsize) {
      $toosmall++;
      next;
    }
    if ($full) {
      print $seq->id, "\t", $seq->alphabet, "\t", $L;
      print "\t", $seq->desc if $desc;
    }
    
    # count how many Ns and -s in all these sequences
    my $s = $seq->seq;
    my $n=0;
    if (not $seq->length) {
      printf STDERR "Warning: zero length sequence '%s'\n", $seq->id;
    }
    else {
      $n = $s =~ s/N/N/gi;
      $nn += $n;
      $n = $s =~ s/-/-/gi;
      $xx += $n;
    }
    
    if ($bases and $seq->alphabet eq 'dna' and $full) {
      my $s = $seq->seq;
      for my $x (qw(N A T G C)) {
        my $n = $s =~ s/$x/$x/gi;
        $n ||= 0;
        printf "\t|%s %6d %4.1f%%", $x, $n, ($n*100/$L);
      }
    }
    print "\n" if $full;
    $numseq++;
    $avglen+=$L;
    $maxlen = $L if $L > $maxlen;
    $minlen = $L if $L < $minlen;
    push @len, $L; # for N50
  }

  @len = sort { $a <=> $b } @len;
  my $cum = 0;
  my $n50 = 0;
  for my $i (0 .. $#len) {
    $cum += $len[$i];
    if ($cum >= $avglen/2) {
  #    $n50 = int( ($len[$i] + $len[$i-1])/2 );
      $n50 = $len[$i];
      last;
    }
  }

  if ($numseq > 0) {
    if($tabbed){
		printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", 
			$fname, $numseq, $avglen, $avglen-$nn-$xx, $nn, $xx, $minlen, $avglen/$numseq, $maxlen, $n50;
	}
	else {
                my $fname2 = shorten_name($fname, $width);
		printf  "%-${width}.${width}s no=%d bp=%d ok=%d Ns=%d gaps=%d min=%d avg=%d max=%d N50=%d", 
			$fname2, $numseq, $avglen, $avglen-$nn-$xx, $nn, $xx, $minlen, $avglen/$numseq, $maxlen, $n50;
	}
    print " (minsize=$minsize, skipped $toosmall)" if $minsize > 0;
    print "\n";
  }
}

#----------------------------------------------------------------------

sub shorten_name {
  my($s, $w) = @_;
  return $s if $w >= length($s);
#  $s =~ s/\..*$//;
  if (length($s) >= $w) {
    $s =~ s/[aeiou]//gi;
  }
  $s = substr $s, -$w;
  return $s;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"desc!",  VAR=>\$desc, DEFAULT=>0, DESC=>"Description at end"},
    {OPT=>"bases!",  VAR=>\$bases, DEFAULT=>0, DESC=>"Base composition (implies --full)"},
    {OPT=>"full!",  VAR=>\$full, DEFAULT=>0, DESC=>"Print details for each sequence"},
    {OPT=>"each!",  VAR=>\$each, DEFAULT=>0, DESC=>"Don't combine results of all input files into one"},
    {OPT=>"minsize=i",  VAR=>\$minsize, DEFAULT=>0, DESC=>"Minimum size to include in calcs"},
    {OPT=>"w|width=i",  VAR=>\$width, DEFAULT=>25, DESC=>"Max. width of filename column"},
    {OPT=>"t|tabbed!",  VAR=>\$tabbed, DEFAULT=>0, DESC=>"Produce tab delimited output table"},
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
  print "Usage: $0 [options] file.fasta [ more.fasta ... ]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
