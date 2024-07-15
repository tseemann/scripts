#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Tools::Run::Alignment::MAFFT;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::TCoffee;
use Bio::Tools::Run::Alignment::Probcons;

my(@Options, $verbose, $program, $format);
setOptions();

my %seq;
my $in = Bio::SeqIO->new(-file=>$ARGV[0], -format=>'Fasta');
while (my $seq = $in->next_seq) {
  my $id = $seq->display_id;
  die "already seen sequence '$id' in $ARGV[0]" if exists $seq{$id};
  $seq{$id} = $seq;
} 

my @param = (-diags=>undef, -stable=>undef);
my $factory = "Bio::Tools::Run::Alignment::$program"->new(@param);
$factory->quiet(1) if $factory->can('quiet');
$factory->stable(1) if $factory->can('stable');

my $bls = Bio::SearchIO->new(-file=>$ARGV[1]); #, -format=>'blast');
while (my $res = $bls->next_result) {
  print STDERR "Processing: ", $res->query_name,"\n";  
  print "\n",$res->query_name,"\n\n";  
  my @hit;
  while (my $hit = $res->next_hit) {
#    print "\t", $hit->name,"\n";
    push @hit, $seq{ $hit->name };
    print $hit[-1]->display_id, " ", $hit[-1]->desc, "\n";
  }
  if (@hit) {
    @hit = sort { $a->display_id cmp $b->display_id } @hit;
    my $align = $factory->align(\@hit);
    my $aout = Bio::AlignIO->new(-fh=>\*STDOUT,-format=>$format);
    $aout->write_aln($align);
  }
  else {
    print "NO HITS\n";
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"p|program=s",  VAR=>\$program, DEFAULT=>'MAFFT', DESC=>"Aligner: Muscle MAFFT Clustalw TCoffee Probcons"},
    {OPT=>"f|format=s",  VAR=>\$format, DEFAULT=>'clustalw', DESC=>"Output alignment format: clustalw stockholm msf fasta nexus"},
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
  print "Usage: $0 [options] <database.fasta> <report.blast>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
