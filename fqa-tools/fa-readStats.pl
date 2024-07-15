#!/usr/bin/env perl
#       readStats.pl
#       
#       Copyright 2011 Simon Gladman <gla048@boomer-wb>
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
use warnings;
#use FindBin;
#use lib "$FindBin::Bin";
use Data::Dumper;
use IO::File;

my(@Options, $verbose, $oneline, $each, $minlen, $header);
setOptions();

my %inputs = $each ? (map { ($_ => IO::File->new($_)) } @ARGV) : ( '(stdin)' => \*ARGV ) ;

if($oneline && $header){
	print "Filename\tyield\tmin\tmedian\tmax\taverage\tn50\n";
}

my $THROBBER = 234655;

for my $fname (sort keys %inputs) {
	
	#print STDERR "$fname\n";
	my $sum = 0;
	my $count = 0;
	
	my $max = 0;
	my $min = 1E100;
	my %lengths;
	my $l = 0;
	my $io = $inputs{$fname};
	my $calc_count = 0;
	
	while(<$io>){
	#while(<ARGV>){
		chomp;
		$count ++;
		if(/^>/){
			next if $count == 1;
			if($l > $minlen){
				$calc_count ++;
				
				$sum += $l;
				$max = $l if $l > $max;
				$min = $l if $l < $min;
				$lengths{$l} ++;
			}
			$l = 0;
			
		}
		else {
			$l += length($_);
			
		}
		if($count % $THROBBER == 0){
			print STDERR "\r$count reads processed" if $verbose; 
		}
	}
	print STDERR "\n" if $verbose;
	
	my $avg = $sum/$calc_count;
	
	my $med_count = 0;
	my $med;
	foreach my $l (sort {$a <=> $b} keys %lengths){
		$med_count += ($lengths{$l}) || 0;
		if($med_count >= $calc_count/2){
			$med = $l;
			last;
		}
	}
	
	my $n50_count = 0;
	my $n50;
	foreach my $l (sort {$a <=> $b} keys %lengths){
		$n50_count += ($l * $lengths{$l}) || 0;
		if($n50_count > $sum/2){
			$n50 = $l;
			last;
		}
	}
	
	
	if($oneline){
		printf "$fname\t$sum\t$calc_count\t$min\t$med\t$max\t%.2f\t$n50\n", $avg;
	}
	else{
		print "\nFile:\t$fname\nNumber of reads:\t$calc_count\n";
		print "Yield:\t$sum\n";
		printf "Average length:\t%.2f\n", $avg;
		print "Minimum length:\t$min\n";
		print "Maximum length:\t$max\n";
		print "Median length:\t$med\n";
		print "n50 length:\t$n50\n";
	}
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"o|oneline!", VAR=>\$oneline, DEFAULT=>0, DESC=>"One line output"},
		{OPT=>"e|each!", VAR=>\$each, DEFAULT=>0, DESC=>"Do each file individually"},
		{OPT=>"m|minlen=i", VAR=>\$minlen, DEFAULT=>0, DESC=>"The minimum length contig/read to include in calcs."},
		{OPT=>"t|tabbed_header!", VAR=>\$header, DEFAULT=>0, DESC=>"Print a header above the one line output"},
	);

	#(@ARGV < 1) && (usage());

	&GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

	# Now setup default values.
	foreach (@Options) {
		if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
		${$_->{VAR}} = $_->{DEFAULT};
		}
	}
}

sub usage {
	print "Usage: $0 [options] reads.fa\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
