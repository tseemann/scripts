#!/usr/bin/env perl
use strict;
use Text::CSV;

my(@Options, $verbose, $noheader, $notitle, $sepchar);
setOptions();

# $status = $csv->combine(@columns);    # combine columns into a string
#        $line   = $csv->string();             # get the combined string
#        $status       = $csv->status ();      # get the most recent status
#        $bad_argument = $csv->error_input (); # get the most recent bad argument
#        $diag         = $csv->error_diag ();  # if an error occured, explains WHY

my $csv = Text::CSV->new( { sep_char => $sepchar } );

print "<TABLE BORDER=1 STYLE='border-collapse: collapse; font-family: sans-serif; font-size: smaller; text-align: center;'>\n";
my $count=0;

while (<ARGV>) {
  chomp;
  my $status = $csv->parse($_);
  $status or die $csv->error_diag."\n".$csv->error_input;
  $count++;
  my $TD = ($count==1 and !$noheader) ? 'TH' : 'TD';
  print "<TR>\n";
  my @f = map { defined $_ ? $_ :  '&nbsp;' } $csv->fields;
  print map { "<$TD>$_\n" } @f;
}
print "</TABLE>\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"noheader!",  VAR=>\$noheader, DEFAULT=>0, DESC=>"Don't treat first line as header"},
    {OPT=>"sepchar=s",  VAR=>\$sepchar, DEFAULT=>',', DESC=>"Separator character"},
#    {OPT=>"notitle!",  VAR=>\$notitle, DEFAULT=>0, DESC=>"Don't create a title"},
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
