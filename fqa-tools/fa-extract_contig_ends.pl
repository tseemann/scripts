#!/usr/bin/env perl
#
#	fa-extract_contig_ends.pl
#
#	Simon Gladman - 2009
#

use strict;

use Bio::SeqIO;
use Bio::Seq;

my(@Options, $b, $t, $f);
setOptions();

my $sio = Bio::SeqIO->new(-fh => \*ARGV, -format => 'Fasta');

while(my $s = $sio->next_seq()){
	my $sepn = length($s->seq()) - (2 * $b + 2 * $t);
	my $lid = ">" . $s->id . "_left_$sepn";
	my $lseq = substr($s->seq(), $t, $b);
	my $rid = ">" . $s->id . "_right_$sepn";
	my $rseq = substr($s->seq(), (length($s->seq()) - $t - $b), $b);
	if($f){
		my $sio = Bio::SeqIO->new(-fh => \*STDOUT , -format => 'Fasta');
		my $ls = Bio::Seq->new(-id => $lid, -seq => $lseq);
		my $rs = Bio::Seq->new(-id => $rid, -seq => $rseq);
		$sio->write_seq($ls);
		$sio->write_seq($rs);
	}
	else {
		print "$lid\n";
		print "$lseq\n";
		print "$rid\n";
		print "$rseq\n";
	}
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",	VAR=>\&usage,           DESC=>"This help"},
    {OPT=>"bases=i",VAR=>\$b,	DEFAULT=>50,DESC=>"Number of bases to include in the contig ends"},
    {OPT=>"trim=i",	VAR=>\$t,	DEFAULT=>5,	DESC=>"Number of bases to trim off the ends first"},
	{OPT=>"fasta!",	VAR=>\$f,	DEFAULT=>0, DESC=>"60 bases per line in output or single line per sequence"}
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
  print "Usage: $0 [options] <contigs fasta>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
