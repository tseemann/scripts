#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SearchIO;

my(@Options, $verbose, $nonohits, $nohitstring, $desclen);
setOptions();

my $bls = Bio::SearchIO->new(-fh=>\*ARGV, -format=>'blast');
while (my $res = $bls->next_result) {
  if ($res->no_hits_found or $res->hits <= 0) {
    next if $nonohits;
    print $res->query_name, "\t", $nohitstring, "\n";
  }
  else {
    my $hit = $res->next_hit;
    my $desc = $hit->description;
    $desc = substr $desc, 0, $desclen if $desclen;
    print $res->query_name, "\t", $desc, "\n";
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"n|nonohits!",  VAR=>\$nonohits, DEFAULT=>0, DESC=>"Don't show 'no hits'"},
    {OPT=>"s|nohitstring!",  VAR=>\$nohitstring, DEFAULT=>'', DESC=>"String for 'no hit'"},
    {OPT=>"l|desclen=i",  VAR=>\$desclen, DEFAULT=>60, DESC=>"Chop hit desc. to this length"},
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
