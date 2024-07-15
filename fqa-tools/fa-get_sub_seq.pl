#!/usr/bin/env perl
#       fa-get_sub_seq.pl
#       
#       Copyright 2010 Simon Gladman <gla048@boomer-wb>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.


use strict;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";
use Data::Dumper;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;

my(@Options, $verbose, $start, $end, $name);
setOptions();


#----------------------------------------------------------------------
# Read and Write a Fasta File

my $in   = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
while (my $seq = $in->next_seq) {
	my $outseq;
	if($name){
		next unless $seq->id =~ /$name/;
	}
	if($end > $seq->length){ die "End requested ($end) is outside of sequence length (".$seq->length.")";}
	my $outstr = $seq->subseq($start,$end);
	my $outid = $seq->id() . "_$start" . "_$end";
	$outseq = Bio::Seq->new(-id => $outid, -seq => $outstr);
	$out->write_seq($outseq);	
} 




#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"s|start=i", VAR=>\$start, DEFAULT=>0, DESC=>"The start of the subsection to pull out"},
		{OPT=>"e|end=i", VAR=>\$end, DEFAULT=>10, DESC=>"The end of the subsection to pull out"},
		{OPT=>"n|name=s", VAR=>\$name, DEFAULT=>"", DESC=>"A sequence id to match for"},
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
	print "Usage: $0 [options] input_fasta.fa > output.fa\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
