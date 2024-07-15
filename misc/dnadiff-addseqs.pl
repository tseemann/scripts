#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use List::Util qw(max min);
use Fatal;

my(@Options, $verbose, $prefix, $minlen);
setOptions();

print STDERR "Opening: $prefix.delta\n";
open my $DELTA, '<', "$prefix.delta";
my $line = <$DELTA>;
chomp $line;
my($ref,$qry) = split ' ', $line;
print STDERR "Ref: $ref\nQry: $qry\n";

-r $ref or die "$ref has been moved or deleted!";
-r $qry or die "$qry has been moved or deleted!";
my $out;

$out = Bio::SeqIO->new(-file=>">$prefix.ronly", -format=>'fasta');
process_seq($out, $ref, "$prefix.rdiff", "$prefix.unref");
print STDERR "Wrote: $prefix.ronly\n";

$out = Bio::SeqIO->new(-file=>">$prefix.qonly", -format=>'fasta');
process_seq($out, $qry, "$prefix.qdiff", "$prefix.unqry");
print STDERR "Wrote: $prefix.qonly\n";

sub process_seq {
  my($outfh, $sfname, $dfname, $unfname) = @_;

  my $seq_of = fasta_to_hash($sfname);

  print STDERR "Opening: $dfname\n";
  open my $dfh, '<', $dfname;
  while (<$dfh>) {
    next unless m/^(\S+) \s+ (GAP|BRK) \s+ (\d+) \s+ (\d+) \s+ (\d+) /x; # ignore -ve col 5
    next unless $5 >= $minlen;
    my $seq = $seq_of->{$1}->trunc( min($3,$4), max($3,$4) );
    next if $seq->seq =~ m/^N+$/;
    $seq->id( $seq->id . ".$dfname.$3-$4.${5}bp" );
    $outfh->write_seq( $seq );
  }
  
  print STDERR "Opening: $unfname\n";
  open my $ufh, '<', $unfname;
  while (<$ufh>) {
    next unless m/^(\S+)\s+(\d+)/;
    next unless $2 >= $minlen;
    $outfh->write_seq( $seq_of->{$1} );
  }
}

sub fasta_to_hash {
  my($sfname) = @_;
  my %hash;
  print STDERR "Opening: $sfname\n";
  my $sin = Bio::SeqIO->new(-file=>$sfname, -format=>'Fasta');
  while (my $seq = $sin->next_seq) {
    $hash{$seq->id} = $seq;
  }
  return \%hash;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"p|prefix=s",  VAR=>\$prefix, DEFAULT=>'out', DESC=>"Prefix used by 'dnadiff'"},
    {OPT=>"m|minlen=i",  VAR=>\$minlen, DEFAULT=>0, DESC=>"Minimum different length to output"},
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

#    velvet-estimate-exp_cov.pl
#
#    Estimates the expected k-mer coverage parameter (-exp_cov) for velvetg
#    by finding the mode of coverage distribution as presented in stats.txt
#
#    Copyright (C) 2009 Torsten Seemann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


