#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::SearchIO;

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



my(@Options, $verbose);
setOptions();

my $in   = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
while (my $seq = $in->next_seq) {
  $out->write_seq($seq);
} 

my $gbk = Bio::SeqIO->new(-file=>'input.gbk', -format=>'genbank');
while (my $seq = $gbk->next_seq) {
  for my $f ($seq->get_SeqFeatures) {
    next if $f->primary_tag eq 'source';
    print $f->primary_tag;
#    my($tag) = $f->has_tag('locus_tag') ? 
    my($tag) = $f->get_tag_values('locus_tag');      
    $tag ||= '(no tag)';
    print "\t", $tag; 
    print "\t", $f->location->gff_string;
    print "\n";
  }
}

my $bls = Bio::SearchIO->new(-file=>'out.blast', -format=>'blast');
while (my $res = $bls->next_result) {
  next if $res->no_hits_found;
  print STDERR $res->query_name,"\n";
  while (my $hit = $res->next_hit) {
    print STDERR "\t", $hit->name,"\n";
    while (my $hsp = $hit->next_hsp) {
      print STDERR "\t\t", $hsp->significance,"\n";
    }
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
