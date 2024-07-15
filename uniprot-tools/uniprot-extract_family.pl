#!/usr/bin/env perl
use strict;
use SWISS::Entry;
use SWISS::KW;
#use SWISS::OC;
use Data::Dumper;

my(@Options, $verbose, $term);
setOptions();

# Read an entire record at a time
local $/ = "\n//\n";

my $in=0;
my $out=0;

while (<ARGV>) 
{
  # Read the entry
  my $entry = SWISS::Entry->fromText($_);
  $in++;

  # Only specified organism class
  my $tax = $entry->OCs->list;
  next unless $tax;
  next unless grep { $_ eq $term } @$tax ;

  print STDERR join(';',@$tax),"\n" if $verbose;

  # Print the primary accession number of each entry.
#  print $entry->AC, "\n";
 
  if (0) { 
    # grep ^PE uniprot_sprot.dat  | sort | uniq -c | sort -n
    #  65866 PE   1: Evidence at protein level;
    #  66437 PE   2: Evidence at transcript level;
    # 348167 PE   3: Inferred from homology;
    #  13928 PE   4: Predicted;
    #   1482 PE   5: Uncertain;
    $entry->PE->text =~ m/^(\d+)/;
    my $pe = $1;
    next unless $pe <= 3; # 1=prot_evid,2=EST_evid,3=homolo_evid
  }

  my $prod = 'unknown transcript';
  for my $de ($entry->DEs->elements) {
    if ($de->type eq 'Full') {
      $prod = $de->text;
    }
  }
  
  print ">", $entry->AC, " ", $prod, "\n",
#        $tax->[-1], "\n", 
        $entry->SQs->seq, "\n";
  $out++;
  print STDERR "\rFound $out from $in so far...";

}

print STDERR "\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"t|term=s",   VAR=>\$term, DEFAULT=>'Pasteurellales', DESC=>"Only allow these OC types"},
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
  print "Usage: $0 [options] <uniprot.dat>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
