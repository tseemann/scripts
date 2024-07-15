#!/usr/bin/perl -w
use strict;

my(@Options, $inverse);
setOptions();

#0 readname contigname strand 
#3 contigstart contigend 
#5 readstart readend readlength 
#8 score editstring

print STDERR "Assuming input is sorted by readname... sort -k 1\n";
my $nread=0;
my $nwrote=0;
my %seen;

while (my $line = <ARGV>) {
  chomp $line;
  $nread++;
  print STDERR "\rProcessed: $nread" if $nread%9751==0;
  my @f = split m/\t/, $line;
  next unless @f >= 10;
  next if $seen{ $f[0] }++;
  print $line,"\n";
  $nwrote++;
}
print STDERR "\rRead $nread, wrote $nwrote unambiguous mapped reads.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
#    {OPT=>"v!",  VAR=>\$inverse, DEFAULT=>0, DESC=>"Return only NON-EXACT hits instead"},
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
