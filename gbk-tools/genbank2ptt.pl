#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

# This script takes a GenBank file as input, and produces a
# NCBI PTT file (protein table) as output. A PTT file is 
# a line based, tab separated format with fixed column types.
#
# Written by Torsten Seemann
# 18 September 2006
# Updated on Wed 4 June 2014 to support multi genbank

@ARGV or die "Usage: $0 file.gbk > file.ptt";

my $gbk = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'genbank');

while (my $seq = $gbk->next_seq) {
  my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;

  print $seq->description, " - 0..",$seq->length,"\n";
  print scalar(@cds)," proteins\n";
  print join("\t", qw(Location Strand Length PID Gene Synonym Code COG Product)),"\n";

  for my $f (@cds) {
    my $gi = '-';
    $gi = $1 if tag($f, 'db_xref') =~ m/\bGI:(\d+)\b/;
    my $cog = '-';
    $cog = $1 if tag($f, 'product') =~ m/^(COG\S+)/;
    my @col = (
      $f->start.'..'.$f->end,
      $f->strand >= 0 ? '+' : '-',
      ($f->length/3)-1,
      $gi,
      tag($f, 'gene'),
      tag($f, 'locus_tag'),
      $cog,
      tag($f, 'product'),
    );
    print join("\t", @col), "\n";
  }
}

sub tag {
  my($f, $tag) = @_;
  return '-' unless $f->has_tag($tag);
  return join(' ', $f->get_tag_values($tag));
}

