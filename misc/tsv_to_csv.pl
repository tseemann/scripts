#!/usr/bin/env perl
use strict;
use Text::CSV;

my(@Options, $verbose, $blank);
setOptions();

# $status = $csv->combine(@columns);    # combine columns into a string
#        $line   = $csv->string();             # get the combined string
#        $status       = $csv->status ();      # get the most recent status
#        $bad_argument = $csv->error_input (); # get the most recent bad argument
#        $diag         = $csv->error_diag ();  # if an error occured, explains WHY

my $csv = Text::CSV->new();

my $N = 0;

while (<>) {
  chomp;
  my @col = split m/\t/;
  if ($N==0) {
    $N = @col; # first row;
  }
  elsif ($N != @col) {
    while (@col < $N) { push @col, ''; }
  }
  if ($blank) { @col = map { $_ || $blank } @col; }
  my $status = $csv->combine(@col);  
  $status or die $csv->error_diag."\n".$csv->error_input;
  my $line = $csv->string;
  print $line,"\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"blank=s",  VAR=>\$blank, DEFAULT=>"unknown", DESC=>"Blank entry"},
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
