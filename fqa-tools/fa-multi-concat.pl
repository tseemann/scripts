#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use List::MoreUtils qw(all any);

my(@Options, $verbose, $keepstop, $plain);
setOptions();

my $NUM = scalar(@ARGV);
my @in = map { Bio::SeqIO->new(-file=>$_, -format=>'Fasta') } @ARGV;
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $nout = 0;

while (1) {
  my @seq = map { $_->next_seq } @in;
  
  last if @seq == 0;
  if (@seq < $NUM) {
    print STDERR "Aborted: Could not read a sequence from all $NUM files at #$nout in @ARGV\n";
    exit 1;
  }
      
  my $seq = join( '', map { defined $_ ? $_->seq : () } @seq );
  $seq =~ s/\*//g unless $keepstop;
  
  my $desc = $plain 
           ? ''
	   : '['.scalar(@seq).'] '.join( ' ~~~ ', map { $_->desc } @seq );
  
  $out->write_seq(
    Bio::Seq->new( -id => ++$nout, -desc => $desc, -seq => $seq )
  );
} 

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"keepstop!",  VAR=>\$keepstop, DEFAULT=>0, DESC=>"Don't remove stop codon"},
    {OPT=>"plain!",  VAR=>\$plain, DEFAULT=>0, DESC=>"Don't write fasta description"},
#    {OPT=>"minsize=i",  VAR=>\$minsize, DEFAULT=>0, DESC=>"Minimum sequence length"},
#    {OPT=>"maxsize=i",  VAR=>\$maxsize, DEFAULT=>1E10, DESC=>"Maximum sequence length"},
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
  print "Usage: $0 [options] f1.fa f2.fa f3.fa > joined_across.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
