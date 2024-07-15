#!/usr/bin/env perl
use strict;
use Text::CSV;

my(@Options, $verbose);
setOptions();

# $status = $csv->combine(@columns);    # combine columns into a string
#        $line   = $csv->string();             # get the combined string
#        $status       = $csv->status ();      # get the most recent status
#        $bad_argument = $csv->error_input (); # get the most recent bad argument
#        $diag         = $csv->error_diag ();  # if an error occured, explains WHY

my $csv = Text::CSV->new({ binary => 1, allow_loose_quotes=>1, allow_whitespace=>1 });

while (<>) {
  chomp;
  my $status = $csv->parse($_);
  $status or die $csv->error_diag."\n".$csv->error_input;
  print join("\t", $csv->fields),"\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
