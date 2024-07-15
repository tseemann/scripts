#!/usr/bin/env perl

use strict;
use warnings;
use Bio::AlignIO;

my(@Options, $verbose, $infile, $outfile, $informat, $outformat);
setOptions();

my $in  = Bio::AlignIO->new(-file => $infile, -format => $informat);
my $out = Bio::AlignIO->new(-file => ">$outfile", -format => $outformat);

while ( my $aln = $in->next_aln ) { 
  $out->write_aln($aln); 
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"i|infile=s", VAR=>\$infile, DEFAULT=>"", DESC=>"Input file"},
		{OPT=>"o|outfile=s",VAR=>\$outfile, DEFAULT=>"", DESC=>"Output file name"},
		{OPT=>"f|informat=s",VAR=>\$informat, DEFAULT=>"fasta", DESC=>"Input file format"},
		{OPT=>"x|outformat=s",VAR=>\$outformat, DEFAULT=>"stockholm", DESC=>"Output file name"},
	);

	(@ARGV < 1) && (usage());

	&GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

	# Now setup default values.
	foreach (@Options) {
		if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
		${$_->{VAR}} = $_->{DEFAULT};
		}
	}
}

sub usage {
	print STDERR "Usage: $0 [options]\n";
	foreach (@Options) {
		printf STDERR "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
