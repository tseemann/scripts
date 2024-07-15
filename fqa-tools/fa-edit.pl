#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $pos, $insert, $delete, $change);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

# we only operate on the FIRST sequence in the list
my $seq = $in->next_seq;
print STDERR $seq->id, " old length ", $seq->length, "\n";

# substr EXPR,OFFSET,LENGTH,REPLACEMENT

if ($pos > 0 and $pos <= $seq->length+1) {  # need +1 for appending (insert at end)
  my $s = $seq->seq;
  if ($change) {
    substr $s, $pos-1, length($change), $change;
  }
  elsif ($delete) {
    substr $s, $pos-1, $delete, '';
  }
  elsif ($insert) {
    substr $s, $pos-1, 0, $insert;
  }
  else {
    quit("no valid change specified");
  }
  $seq->seq($s);
}
else {
  quit("invalid --pos $pos : out of valid range 1..".$seq->length+1);
}

print STDERR $seq->id, " new length ", $seq->length, "\n";
$out->write_seq($seq);


sub quit {
  print STDERR "Error: @_\n";
  exit;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"p|pos=i",  VAR=>\$pos, DEFAULT=>0, DESC=>"Position to perform edit"},
    {OPT=>"c|change=s",  VAR=>\$change, DEFAULT=>'', DESC=>"Substitute to these letters"},
    {OPT=>"d|delete=i",  VAR=>\$delete, DEFAULT=>0, DESC=>"Delete this many letters"},
    {OPT=>"i|insert=s",  VAR=>\$insert, DEFAULT=>'', DESC=>"Insert these letters"},
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
  print "Usage: $0 [options] needs_editing.fasta > fixed.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
