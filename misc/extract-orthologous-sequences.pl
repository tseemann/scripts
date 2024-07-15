#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Temp qw(tempdir);

my(@Options, $verbose);
setOptions();

my $dir = tempdir( CLEANUP => 1 );

my $qfile = shift @ARGV;
print STDERR "Query is: $qfile\n";

for my $rfile (@ARGV) {
  print STDERR "Searching against: $rfile\n";
  my $cmd = "glsearch36 -L -m 10 -T 8 -n -E "0 0" -b 1 -d 1 $qfile $rfile";
  print STDERR "Command: $cmd\n";
  open SS, "$cmd |";
  close SS;
}


#my $in   = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
#my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
#while (my $seq = $in->next_seq) {
#  $out->write_seq($seq);
#} 


#my $bls = Bio::SearchIO->new(-file=>'out.blast', -format=>'blast');
#while (my $res = $bls->next_result) {
#  next if $res->no_hits_found;
#  print STDERR $res->query_name,"\n";
#  while (my $hit = $res->next_hit) {
#    print STDERR "\t", $hit->name,"\n";
#    while (my $hsp = $hit->next_hsp) {
#      print STDERR "\t\t", $hsp->significance,"\n";
#    }
#  }
#}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
  print STDERR "Usage: $0 [options] <query.fasta> <contigs1.fna> <contigs2.fna ...\n";
  foreach (@Options) {
    printf STDERR "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
