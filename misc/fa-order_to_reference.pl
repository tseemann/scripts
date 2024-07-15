#!/usr/bin/env perl
#       fa-order_to_reference.pl
#       
#       Copyright 2012 Simon Gladman <gla048@hyperion>
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
#use Bio::SeqIO;
#use Bio::Seq;
#use Bio::SearchIO;
use File::Temp qw (tempfile tempdir);

my(@Options, $verbose, $ref, $ctgf, $out);
setOptions();

#check if Mauve is installed or find the Mauve jar file...
my ($mauve) = qx(which Mauve) or die "Mauve not installed on this machine or not in path";
chomp ($mauve);
$mauve .= ".jar";
print STDERR "Mauve installed: $mauve\n" if $verbose;

#my $mauve = "/bio/sw/mauve_2.3.1/Mauve.jar";

#check java is installed
my ($java) = qx(which java) or die "java not installed on this machine or not in path";
print STDERR "Java installed: $java\n" if $verbose;

#check files exist..
unless (-r $ref) { die "Reference file doesn't exist or is not readable\n$!"; }
unless (-r $ctgf) { die "Contigs file doesn't exist or is not readable\n$!"; }

#make a temp dir
my $dir = tempdir( CLEANUP => 1);

my $cmd = "java -Xmx2000m -cp $mauve org.gel.mauve.contigs.ContigOrderer -output $dir -ref '$ref' -draft '$ctgf' 1>$dir/mauveout.txt";

unless ($verbose) { $cmd .= " 2>/dev/null"; }

print STDERR "$cmd\n" if $verbose;

runcmd($cmd);

open RES, "$dir/mauveout.txt" or die "Couldn't open mauve stdout for reading.\n";
my @files;
while(<RES>){
	chomp;
	if(m/^done:\s+(.*)$/){
		push @files, $1;
		#print STDERR $1 . "\n";
	}
}
close RES;

if ($verbose) {
	foreach(@files){
		print STDERR "$_\n";
	}	
}

my $resultfile = $files[-1];
$resultfile =~ s/_contigs.tab$/.fas/;
$resultfile =~ s/\/contigs\//\//;

print STDERR "Output file required = $resultfile\n" if $verbose;

unless (-r $resultfile) {die "Mauve output file not found.\n$!";}

open IN, $resultfile;
open OUT, ">$out";
while(<IN>){
	print OUT $_;
}
close IN;
close OUT;

#runcmd subroutine
sub runcmd {
	print STDERR "Running: @_\n";
	system(@_)==0 or die("Could not run command:", @_);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"r|ref=s", VAR=>\$ref, DEFAULT=>"", DESC=>"The reference fasta file"},
		{OPT=>"c|contigs=s", VAR=>\$ctgf, DEFAULT=>"/dev/stdin", DESC=>"The contigs multifasta file"},
		{OPT=>"o|outfile=s", VAR=>\$out, DEFAULT=>"/dev/stdout", DESC=>"The final output file"},
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
	print "Usage: $0 [options] \n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
