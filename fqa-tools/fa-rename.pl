#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $mapfile);
setOptions();

open my $MAP, "<$mapfile" or die "bad --map";
my %map;
while (<$MAP>) {
  chomp;
  my($old,$new) = split ' ';
  $map{$old} = $new;
}
print STDERR "Map file contained ", scalar(keys %map), " ID mappings\n";

my $in  = Bio::SeqIO->new(-fh=>\*ARGV,   -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
my $count = 0;
my $changed = 0;

while (my $seq = $in->next_seq) {
  print STDERR "Processing: ", $seq->id, "\n" if $verbose;
  if (my $new = $map{$seq->id}) {
    $seq->id($new);
    $changed++;
  }
  $out->write_seq($seq);
  $count++;
} 

print STDERR "Changed $changed/$count IDs\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"mapfile=s",  VAR=>\$mapfile, DEFAULT=>'', DESC=>"OLD_ID to NEW_ID map file"},
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
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
