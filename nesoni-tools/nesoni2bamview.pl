#!/usr/bin/perl
#       nesoni2bamview.pl
#       
#       Copyright 2011 Simon Gladman CSIRO <simon.gladman@csiro.au>
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
use warnings;

my(@Options, $verbose, $width, $refline, $unmasked, $html, $ansi);
setOptions();

my $dir = $ARGV[0];

-r "$dir/report.txt" or die "report.txt doesn't exist. $!\n";
-r "$dir/reference.fa" or die "reference.fa doesn't exist. $!\n";
-r "$dir/alignments_filtered_sorted.bam" or die "alignments_filtered_sorted.bam doesn't exist. $!\n";

my $file = "$dir/report.txt";

my @reportfile = qx(cat $file);
shift @reportfile;

if($html){
	print "<a name='top'><h1>SNP index</h1></a>\n";
	print '<p><OL>' . "\n";
	my $count = 1;
	foreach (@reportfile){
		chomp;
		my @tmp = split /\s+/, $_;
		my $pos = $tmp[1];
		$tmp[1] = "<a href='#$count'>" . $pos  . "</a>";
		my $indexline = join " ", @tmp;
		print "<li>$indexline\n";
		$count ++;
	}
	print '</ol></p>'
}
my $count = 1;
foreach(@reportfile){
	my $cmd = "bam-view.pl -b '$dir/alignments_filtered_sorted.bam' -f '$dir/reference.fa' -w $width -c $refline ";
	$cmd .= "-u " if $unmasked;
	$cmd .= "--ansi " if $ansi;
	chomp;
	print STDERR "Processing:	$_\n" if $verbose;
	print "\n";
	print "<p><a name='$count'><a href='#top'>Home</a><p><hr><h3>SNP Reported at: $_</h3>" if $html;
	print "\n*********************************************************\n\nSNP Reported at: $_\n\n" if !$html;
	my @tmp = split /\s+/, $_;
	my $sequence = $tmp[0];
	my $pos = $tmp[1];
	my $start = $pos - int($width/2);
	$cmd .= "-i \Q$sequence\E -s $start";
	print "<pre>\n" if $html;
	system($cmd) == 0 or die "system call for $cmd failed! $!\n";
	print "</pre><p><a href='#top'>Home</a>\n" if $html;
	$count ++;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"w|width=i", VAR=>\$width, DEFAULT=>100, DESC=>"The width of the ascii pileup display output"},
		{OPT=>"c|consfreq=i", VAR=>\$refline, DEFAULT=>0, DESC=>"Print consensus every N lines (0=just at top)"},
		{OPT=>"u|unmasked!", VAR=>\$unmasked, DEFAULT=>0, DESC=>"Show all bases, don't use '.' for same"},
		{OPT=>"html!", VAR=>\$html, DEFAULT=>0, DESC=>"Make a nice html output."},
		{OPT=>"ansi!", VAR=>\$ansi, DEFAULT=>0, DESC=>"Make a nice ansi output for certain terminals."},
	);

	(@ARGV < 1) && (usage());

	&GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();
	
	die "--ansi and --html are incompatible. please choose only 1\n" if $ansi && $html;
	
	
	# Now setup default values.
	foreach (@Options) {
		if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
		${$_->{VAR}} = $_->{DEFAULT};
		}
	}
}

sub usage {
	print "Usage: $0 [options] nesoni_dir > output.txt\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
