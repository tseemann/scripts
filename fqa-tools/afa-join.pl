#!/usr/bin/env perl
use warnings;
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Align::DNAStatistics;
use Data::Dumper;

my(@Options, $verbose, $informat, $outformat, $out_fn);
setOptions();

my %out;

for my $fname (@ARGV) {
  print STDERR "Open: $fname\n";
  my $in = Bio::AlignIO->new(-file=>$fname, -format=>$informat);
  my $aln = $in->next_aln;
  if (not $aln->is_flush) {
    print STDERR "\tnot flush!\n";
    next;
  }
  print STDERR "\tappending sequences\n";
  for my $seq ($aln->each_seq) {
    $out{ $seq->id } .= $seq->seq;
  }      
}

print STDERR "Sequences: ", scalar(keys %out), "\n";

print STDERR "Building alignment...\n";
my $aln = Bio::SimpleAlign->new;
for my $id (sort keys %out) {
  my $s = $out{$id};
  my $s2 = $s;
  $s2 =~ s/-//g;
  printf STDERR "$id\t%d\t%d\n", length($s), length($s2);
  $aln->add_seq( Bio::LocatableSeq->new(-id=>$id, -seq=>$out{$id}, -start=>1, -end=>length($s2)) );
#  $aln->add_seq( Bio::Seq->new(-id=>$id, -seq=>$out{$id}) );
}

print STDERR "Removing gaps...\n";
$aln->remove_gaps;

print STDERR "Writing alignment: $out_fn\n";
my $out = Bio::AlignIO->new(-file=>">$out_fn", -format=>$outformat);
$aln->set_displayname_flat;
$out->write_aln($aln);

#my $stats = Bio::Align::DNAStatistics->new;
#my $kak = $stats->calc_all_KaKs_pairs($aln);
#print Dumper($kak);

print STDERR "Done.\n";
            
#.......................................................................
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"if=s",  VAR=>\$informat, DEFAULT=>'fasta', DESC=>"Input format"},
    {OPT=>"of=s",  VAR=>\$outformat, DEFAULT=>'fasta', DESC=>"Output format"},
    {OPT=>"out=s",  VAR=>\$out_fn, DEFAULT=>'/dev/stdout', DESC=>"Output alignment file"},
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
