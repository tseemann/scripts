#!/usr/bin/env perl

use strict;
use Bio::Seq;
use Bio::SeqIO;

open IN, $ARGV[0];

my $count = 1;

while(<IN>){
	chomp;
	if(/^>/){
		s/>/>contig_$count /;
		$count ++;
	}
	print $_ . "\n";
}
