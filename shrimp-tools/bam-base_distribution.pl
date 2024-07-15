#!/usr/bin/perl
use strict;

# per-base nucleotide frequences from a BAM file
# torsten@seemann.id.au

my %len_of;

for my $bamfn (@ARGV) {
  -r $bamfn or die "can't read bam file: $bamfn";
  foreach (qx(samtools view -H \Q$bamfn\E)) {
    next unless m/^\@SQ\s+SN:(\S+)\s+LN:(\d+)$/;
    $len_of{$1} = $2;
  }
#  print Dumper(\%len_of);
  
  # print header line
  print join("\t", qw(SEQUENCE POSITION A C G T)),"\n";
  
  # open BAM file
  open PILE, "samtools mpileup -d 1000 \Q$bamfn\E 2> /dev/null |";
  while (<PILE>) {
    # parse a line
    chomp;
    my($id, $pos, undef, $depth, $bases) = split m/\t/;
    # count up occurrences of all chars
    my %freq;
    map { $freq{$_}++ } (split m//, uc($bases));
    # print out frequencies
    print join("\t", 
      $id, 
      $pos, 
      $freq{'A'} || 0, 
      $freq{'C'} || 0, 
      $freq{'G'} || 0, 
      $freq{'T'} || 0,
    ), "\n";
  }              
}

