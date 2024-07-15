#!/home/linuxbrew/.linuxbrew/bin/perl -w
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $format, $moveid);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $count=0;
while (my $seq = $in->next_seq) {
  $count++;
#  print STDERR "\rProcessing: $seq->display_id [$count]";
  if ($moveid) {
    $seq->desc( join ' ', $seq->display_id, $seq->desc );
  }
  $seq->display_id( sprintf $format, $count );
  $out->write_seq($seq);
} 
#print STDERR "\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"format=s",  VAR=>\$format, DEFAULT=>'gnl|VBC|scaffold%04d', DESC=>"Format: %d=seqnum (%04d allowed too etc.)"},
    {OPT=>"moveid!",  VAR=>\$moveid, DEFAULT=>0, DESC=>"Move old ID into Description field"},    
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
  print "Usage: $0 [options] original.fasta > with_new_IDs.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
