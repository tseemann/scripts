#!/usr/bin/env perl
use strict;
use Bio::SearchIO;
use List::Util qw(min max);
use POSIX qw(ceil floor);

my(@Options, $verbose, $format, $width, $namewidth, $pcid, $evalue);
setOptions();

my $idw = $namewidth;
$width ||= 80;

$pcid /= 100;  

my $in = Bio::SearchIO->new(-fh=>\*ARGV, -format=>$format);

while (my $bls = $in->next_result) {
  print STDERR "Showing report...\n" if $verbose;
  show_report($bls);
}

sub show_report {
  my($bls) = @_;
  return if $bls->hits <= 0;
  print ">",$bls->query_name," ",$bls->query_description," ", $bls->query_length, "\n";
  die "query length was zero!" if $bls->query_length <= 0;
  my $chunk = int(ceil($bls->query_length / $width));
  print " "x($idw-1);
  print "1 ";
  print ("."x$width);
#  print " Eval\n";
  my @line;
  print " ",$bls->query_length,"\n";
  while (my $hit = $bls->next_hit) {
    while (my $hsp = $hit->next_hsp) {
      push @line, [ $hit->name, $hsp ];
    }
  }
  @line = sort { $a->[1]->start('qry') <=> $b->[1]->start('qry') } @line;
  
  for my $line (@line) {
    my($hitname, $hsp) = @$line;
    next if $hsp->significance > $evalue;
    next if $hsp->frac_identical < $pcid;
    my $begin = min( $hsp->start('qry'), $hsp->end('qry') );
    my $indent = int(floor( $begin / $chunk));
    my $len = int(ceil($hsp->length / $chunk));
    my $dir = $hsp->strand('hit') > 0 ? '>' : '<';
    printf "%-${idw}.${idw}s ", $hitname;
    print " "x$indent;
    print ${dir}x$len;
    print " "x($width-$indent-$len);
    printf " %3.0f%%", 100*$hsp->frac_identical;
    print " ",$hsp->significance;
    print " ",$hsp->length;
    print "\n";
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"f|format=s",  VAR=>\$format, DEFAULT=>'blast', DESC=>"Format: blast blastxml"},
#    {OPT=>"s|sort!",  VAR=>\$sort, DEFAULT=>0, DESC=>"Sort hits by position"},
    {OPT=>"w|width=i",  VAR=>\$width, DEFAULT=>($ENV{'COLUMNS'} ||  80), DESC=>"Output width for coverage "},
    {OPT=>"n|namewidth=i",  VAR=>\$namewidth, DEFAULT=>25, DESC=>"Output width for name"},
    {OPT=>"e|evalue=f",  VAR=>\$evalue, DEFAULT=>1, DESC=>"Evalue cutoff"},
    {OPT=>"i|pcid=f",  VAR=>\$pcid, DEFAULT=>0, DESC=>"Percent ID cutoff"},
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
  print "Usage: $0 [options] < blast.output\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
