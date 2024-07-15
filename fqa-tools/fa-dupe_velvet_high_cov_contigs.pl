#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use List::Util qw(max);

my(@Options, $verbose, $minlen, $mincov, $maxlen);
setOptions();

my %s;

my $in  = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
while (my $seq = $in->next_seq) {
  # >NODE_42_length_272963_cov_54.489731
  next unless $seq->id =~ m/^NODE_(\d+)_length_(\d+)_cov_([\d\.]+)/;
  next if $2 < $minlen;
  next if $3 < $mincov;
  $s{$1} = {LEN => $2, COV => $3, SEQ => $seq };
  printf STDERR "Read: %s\n", $seq->id if $verbose;
} 
printf STDERR "Accepted %d contigs.\n", scalar(keys %s);

my $top=0;
my $bot=0;
for my $id (keys %s) {  
  $top += $s{$id}{LEN} * $s{$id}{COV};
  $bot += $s{$id}{LEN};
}
my $acov = $top / $bot;
printf STDERR "Estimated mean coverage: %.2f\n", $acov;

my $nout=0;

my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
#for my $id (sort { $s{$a}{COV} <=> $s{$b}{COV} } keys %s) {  
for my $id (sort { $s{$b}{LEN} <=> $s{$a}{LEN} } keys %s) {  
  my $seq = $s{$id}{SEQ};
  my $copies = int( $s{$id}{COV} / $acov ) || 1; # err on conservative side
  print STDERR "[$id] $s{$id}{COV} / $acov => $copies copies\n" if $verbose;
  if ($copies > 1 and $seq->length <= $maxlen) {
    for my $c (1 .. $copies) {
      $seq->id( sprintf 'NODE_%s_length_%d_cov_%.2f_copy_%d', 
        $id, $s{$id}{LEN}, $s{$id}{COV}/$copies, $c );
      $out->write_seq($seq);
      $nout++;
    }
  }
  else {
    $out->write_seq($seq);
    $nout++;
  }
}
print STDERR "Wrote $nout contigs.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"c|mincov=f",  VAR=>\$mincov, DEFAULT=>3, DESC=>"Minimum coverage to allow"},
    {OPT=>"l|minlen=f",  VAR=>\$minlen, DEFAULT=>200, DESC=>"Minimum length to allow"},
    {OPT=>"m|maxlen=i",  VAR=>\$maxlen, DEFAULT=>20000, DESC=>"Maximum length contig to duplicate (beware of plasmids!)"},
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
  print "Usage: $0 [options] <contigs.fa>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
