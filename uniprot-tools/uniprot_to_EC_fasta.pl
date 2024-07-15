#!/usr/bin/env perl
#       swiss-get_ECs.pl
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

#	Modified by Torsten Seemann 1 July 2010
#	- use STDIN (so .gz works)
#	- changed ID to sprot_id (so unique)
#	- removed Bio::SeqIO dependency

use strict;
use Getopt::Long;
use Data::Dumper;

my(@Options, $verbose);
setOptions();

#----------------------------------------------------------------------
# Read and Write a Fasta File

my $x;
my $count = 0;
my $record_count = 0;

# need to replace this loop with $/ = '//' separator

while(<>){
	if( m{^//} ){
		process_record($x);
		$x = "";
		$record_count++;
	}
	else{
		$x .= $_;
	}
}

print STDERR "$count records with EC numbers found..\n";
print STDERR "$record_count total records\n";



sub process_record {
	my $r = shift;
	if($r =~ m/DE\s+EC=([\d-]+\.[\d-]+\.[\d-]+\.[\d-]+)/){
		my $EC = $1;
		my $desc = "";
		$count ++;
		#now get the DE RecName Full=
		if($r =~ m/DE\s+RecName:\s+Full=(.*);/){
			$desc = $1;
		}
		my $seq = "";
		my @temp = split "\n", $r;
		foreach my $line (@temp) {
			if($line =~ m/^\s+/){
				$seq .= $line;
			}
		}
		$seq =~ s/\s+//g;
		$r =~ m/ID   (\w+)/;
		my $ID = $1 || "noid_$count";
		print ">$ID $EC $desc\n$seq\n";
#		my $bioseq = Bio::Seq->new(-id => $EC, -desc => $desc, -seq => $seq);
#		$out->write_seq($bioseq);
	}
	return;	
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
	);

#	(@ARGV < 1) && (usage());

	&GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

	# Now setup default values.
	foreach (@Options) {
		if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
		${$_->{VAR}} = $_->{DEFAULT};
		}
	}
}

sub usage {
	print "Usage: $0 [options] < uniprot_sprot.dat > Enzymes.faa\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
