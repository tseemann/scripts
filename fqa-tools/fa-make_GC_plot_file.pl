#!/usr/bin/env perl
#       fa-make_GC_plot_file.pl
#       
#       Copyright 2011 Simon Gladman <simon.gladman@csiro.au>
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
use Bio::SeqIO;
use Bio::Seq;
use List::Util qw(sum);

my(@Options, $verbose, $window, $delimeter, $spars);
setOptions();

my $THROBBER = 56789;

#----------------------------------------------------------------------
# Read and Write a Fasta File

my $in   = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');

while (my $seq = $in->next_seq) {
	my $s = $seq->seq();
	print STDERR "Working on " . $seq->id . " with:\nWindow size: $window\nSparsity: $spars\n" if $verbose;
	
	my $pos = 0;
	my $maxpos = length($s) - 1;
	my $halfwindow = int($window/2);
	my $maxfullwindow = $maxpos - $halfwindow;
	my $maxgc = 0;
	my $mingc = 100;
	my $gcsum = 0;
	my $sum = 0;
	my $length;
	my $gc;
	print STDERR "Sequence length: $maxpos\nTranslating\n" if $verbose;
	
	#convert the entire sequence to 0's and 1's.  (0 for A&T, 1 for G&C)
	$s =~ tr/aAtTcCgG/00001111/;
	my @tmp = split //, $s;
	
	print STDERR "Finished translating\n" if $verbose;
	
	$sum = sum(@tmp[0..$halfwindow-1]);
	$length = $halfwindow;
	$gc = $sum/$length *100;
	$maxgc = $gc;
	$mingc = $gc;
	$gcsum = $gc;
	
	my @outcoords;
	my @outgcs;
	
	printf ("%d$delimeter%.1f\n",0,$gc);
	
	for(my $i = 1; $i < $halfwindow; $i ++){
		$length = $i + $halfwindow;
		$sum += $tmp[$i+$halfwindow];
		$gc = $sum/$length * 100;
		$maxgc = $gc if $gc > $maxgc;
		$mingc = $gc if $gc < $mingc;
		$gcsum += $gc;
		if($i % $spars == 0){
			push @outcoords, $i;
			push @outgcs, $gc;
		}
	}
	
	
	for(my $i = $halfwindow; $i <= $maxfullwindow; $i ++){
		if($i % $THROBBER == 0){
			print STDERR "\r$i" if $verbose;
		}		
		$length = $window;
		$sum = $sum + $tmp[$i+$halfwindow] - $tmp[$i - $halfwindow];
		$gc = $sum/$length * 100;
		$maxgc = $gc if $gc > $maxgc;
		$mingc = $gc if $gc < $mingc;
		$gcsum += $gc;
		if($i % $spars == 0){
			push @outcoords, $i;
			push @outgcs, $gc;
		}
	}
	
	for(my $i = $maxfullwindow + 1; $i <= @tmp; $i ++){
		$length = $maxpos - $i + $halfwindow;
		$sum -= $tmp[$i - $halfwindow];
		$gc = $sum/$length * 100;
		$maxgc = $gc if $gc > $maxgc;
		$mingc = $gc if $gc < $mingc;
		$gcsum += $gc;
		if($i % $spars == 0){
			push @outcoords, $i;
			push @outgcs, $gc;
		}
	}
	
	print STDERR "\nWriting out file\n" if $verbose;
	
	my $avggc = $gcsum/$maxpos;
	
	printf STDERR ("Max GC = %.1f\n", $maxgc) if $verbose;
	printf STDERR ("Min GC = %.1f\n",$mingc) if $verbose;
	printf STDERR ("Avg GC = %.1f\n", $avggc) if $verbose;
	
	for(my $i = 0; $i < @outcoords; $i ++){
		printf ("%d$delimeter%.1f\n",$outcoords[$i],$outgcs[$i]); 
	}
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"w|window=i", VAR=>\$window, DEFAULT=>500, DESC=>"The window size to calculate GC content from"},
		{OPT=>"s|sparsity=i", VAR=>\$spars, DEFAULT=>10, DESC=>"One point in the plot every X bases..."},
		{OPT=>"d|delimeter=s", VAR=>\$delimeter, DEFAULT=>"\t", DESC=>"The delimeter to use between the position and the gc content (, or \\t usually)"},
		
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
	print "Usage: $0 [options] input.fa > gcontent.plotfile\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
