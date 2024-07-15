#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Bio::SeqUtils;

my(@Options, $debug, $informat, $outformat);
setOptions();

# read in all files

my @seq;
for my $file (@ARGV) {
  print STDERR "Loading: $file\n";
  my $seqi = Bio::SeqIO->new(-file=>$file, -format=>$informat);
  while (my $seq = $seqi->next_seq) {
    print STDERR "\t",$seq->id,"\n";
    push @seq, $seq;
  }
}

if (not @seq) {
  print STDERR "No sequences to join!\n";
  exit;
}

# concat files

print STDERR "Joining ",scalar(@seq)," sequences...\n";
my $seq = shift(@seq);
Bio::SeqUtils->cat($seq, @seq) or die "could not join sequences";

# write it out

print STDERR "Writing joined sequence to STDOUT...\n";
my $seqo = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>$outformat);
$seqo->write_seq($seq);

print STDERR "Done.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"if|informat=s",  VAR=>\$informat, DEFAULT=>'genbank', DESC=>"Output format"},
    {OPT=>"of|outformat=s",  VAR=>\$outformat, DEFAULT=>'genbank', DESC=>"Output format"},
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
  print "Usage: $0 [options] file1.gb multifile2.gb ... > joined.gb\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
