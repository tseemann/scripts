#! /usr/bin/perl

# Makes a tab-separated file from 
# the 454NewblerMetrics.txt file
# from a newbler assembly
# tested on newbler v 2.3 and 2.5.3
# on both shotgun, shotgun + paired end and transcriptome assemblies
# version 1, May 2011
# by Lex Nederbragt, lex.nederbragt@bio.uio.no

# run as
# newblermetrics.pl 454Newblermetrics.txt
# newblermetrics.pl /path/to/454Newblermetrics.txt

# or
# perl newblermetrics.pl 454Newblermetrics.txt
# perl newblermetrics.pl /path/to/454Newblermetrics.txt


use strict;
use warnings;

###############################
# variables
###############################

my $metrics;		# holds the entire 454NewblerMetrics.txt file
my $section = "";	# section of the file, e.g. rundata
my $level2;			# all lines with a single tab
my $level3;			# all lines with two tabs
my %metrics = ();	# hash with extracted results
my @lib_names;		# names for paired end libraries

###############################
# test inputfile
###############################

# file given?
if (!$ARGV[0]){
	print STDERR "Please add a 454Newblermetrics.txt file on the command line...\n";
	exit[0];
}

# file exists and is a file?
unless (-e $ARGV[0] && -f $ARGV[0]){
	print STDERR "File '$ARGV[0]' does not exist or is not a file...\n";
	exit[0];
} 

# file can be opened?
open METRICS , "<$ARGV[0]" or die "File '$ARGV[0]' can't be opened:\n$!";


# read in the file
$/=undef; # set the record to 'slurp' the file 
$metrics = <METRICS>;

# correct file type?
unless ($metrics =~ /454 Life Sciences Corporation/ && $metrics =~ /Newbler Metrics Results/ ){
	print STDERR "File '$ARGV[0]' does not appear to be a 454NewblerMetrics file...\n";
	exit[0];
}

if ($metrics =~ /Date of Mapping: /){
	print STDERR "The script currently only works on 454NewblerMetrics.txt files
from newbler assemblies,
not from mappings (gsMapper, runmapping)...\n";
	exit[0];
}


###############################
# process inputfile
###############################

foreach (split /\n/, $metrics){
	$section = $_ if /^\w/;

	# runData/pairedReadData
	if ($section eq "runData" || $section eq "pairedReadData"){
		if( /(numberOf.+) = (\d+), (\d+);/){
			$metrics{'reads'}{"$1Raw"}+=$2;
			$metrics{'reads'}{"$1Trimmed"}+=$3;
		}
	next;
	}
	
	# consensusResults section
	if ($section eq "consensusResults"){
		# type of metric/status is on level 2
		$level2 = $1 if /^\t(\w+)/;
		
		# pairedReadStatus
		push @lib_names, $1 if/libraryName\s+= "(.+)";/;
		$metrics{$lib_names[-1]}{$1}=$2 if /(pairDistance...)\s+= ([0-9\.]+);/;
		
		# other metrics
		$metrics{$level2}{$1}=$2 if /^\t\t(\w+)\s+= ([0-9\.]+)/;
		next;
	}
}
$/="\n"; # reset the record separator

# fix spelling mistake 'bug' from newbler 2.3
if ($metrics{'isotigMetrics'}{'numberWithOneConitg'}){
	$metrics{'isotigMetrics'}{'numberWithOneContig'} = $metrics{'isotigMetrics'}{'numberWithOneConitg'}
}


###############################
# output
###############################

print "Input\n";
print "Number of reads\t", $metrics{'reads'}{'numberOfReadsRaw'}, "\n";
print "Number of bases\t", $metrics{'reads'}{'numberOfBasesRaw'}, "\n";
print "Number of reads trimmed\t", $metrics{'reads'}{'numberOfReadsTrimmed'}, "\t",
	sprintf ("%.1f",	100*$metrics{'reads'}{'numberOfReadsTrimmed'}/
						$metrics{'reads'}{'numberOfReadsRaw'}), "%\n";
print "Number of bases trimmed\t", $metrics{'reads'}{'numberOfBasesTrimmed'}, "\t",
	sprintf ("%.1f",	100*$metrics{'reads'}{'numberOfBasesTrimmed'}/
						$metrics{'reads'}{'numberOfBasesRaw'}), "%\n";
print "\n";

print "Consensus results\n";
print "Number of reads assembled\t", $metrics{'readStatus'}{'numberAssembled'},"\t",
	(sprintf "%.1f",	100*$metrics{'readStatus'}{'numberAssembled'}/
						$metrics{'reads'}{'numberOfReadsTrimmed'})."%\n";
print "Number partial\t", $metrics{'readStatus'}{'numberPartial'},"\t",
	(sprintf "%.1f",	100*$metrics{'readStatus'}{'numberPartial'}/
						$metrics{'reads'}{'numberOfReadsTrimmed'})."%\n";
print "Number singleton\t", $metrics{'readStatus'}{'numberSingleton'},"\t",
	(sprintf "%.1f",	100*$metrics{'readStatus'}{'numberSingleton'}/
						$metrics{'reads'}{'numberOfReadsTrimmed'})."%\n";
print "Number repeat\t", $metrics{'readStatus'}{'numberRepeat'},"\t",
	(sprintf "%.1f",	100*$metrics{'readStatus'}{'numberRepeat'}/
						$metrics{'reads'}{'numberOfReadsTrimmed'})."%\n";
print "Number outlier\t", $metrics{'readStatus'}{'numberOutlier'},"\t",
	(sprintf "%.1f",	100*$metrics{'readStatus'}{'numberOutlier'}/
						$metrics{'reads'}{'numberOfReadsTrimmed'})."%\n";
print "Number too short\t", $metrics{'readStatus'}{'numberTooShort'},"\t",
	(sprintf "%.1f",	100*$metrics{'readStatus'}{'numberTooShort'}/
						$metrics{'reads'}{'numberOfReadsTrimmed'})."%\n";
						
print "\n";


if (exists $metrics{'scaffoldMetrics'}{'numberOfScaffolds'}){
	print "Scaffold Metrics\n";
	print "Number of scaffolds\t", $metrics{'scaffoldMetrics'}{'numberOfScaffolds'}, "\n";
	print "Number of bases\t", $metrics{'scaffoldMetrics'}{'numberOfBases'}, "\n";
	print "Average scaffold size\t", $metrics{'scaffoldMetrics'}{'avgScaffoldSize'}, "\n";
	print "N50 scaffold size\t", $metrics{'scaffoldMetrics'}{'N50ScaffoldSize'}, "\n";
	print "Largest scaffold size\t", $metrics{'scaffoldMetrics'}{'largestScaffoldSize'}, "\n";
	print "\n";
}

if (exists $metrics{'isogroupMetrics'}{'numberOfIsogroups'}){
 	print "Isogroup Metrics\n";
 	print "Number of isogroups\t", $metrics{'isogroupMetrics'}{'numberOfIsogroups'}, "\n";
 	print "Average contig count\t", $metrics{'isogroupMetrics'}{'avgContigCnt'}, "\n";
 	print "Largest contig count\t", $metrics{'isogroupMetrics'}{'largestContigCnt'}, "\n";
 	print "Number with one contig\t", $metrics{'isogroupMetrics'}{'numberWithOneContig'}, "\n\n";
 	print "Average isotig count\t", $metrics{'isogroupMetrics'}{'avgIsotigCnt'}, "\n";
 	print "Largest isotig count\t", $metrics{'isogroupMetrics'}{'largestIsotigCnt'}, "\n";
 	print "Number with one isotig\t", $metrics{'isogroupMetrics'}{'numberWithOneIsotig'}, "\n\n";

 	print "Isotig Metrics\n";
 	print "Number of Isotigs\t", $metrics{'isotigMetrics'}{'numberOfIsotigs'}, "\n";
 	print "Average contig count\t", $metrics{'isotigMetrics'}{'avgContigCnt'}, "\n";
 	print "Largest contig count\t", $metrics{'isotigMetrics'}{'largestContigCnt'}, "\n";
 	print "Number with one contig\t", $metrics{'isotigMetrics'}{'numberWithOneContig'}, "\n\n";
 	print "Number of bases\t", $metrics{'isotigMetrics'}{'numberOfBases'}, "\n";
 	print "Average isotig size\t", $metrics{'isotigMetrics'}{'avgIsotigSize'}, "\n";
 	print "N50 isotig size\t", $metrics{'isotigMetrics'}{'N50IsotigSize'}, "\n";
 	print "Largest isotig\t", $metrics{'isotigMetrics'}{'largestIsotigSize'}, "\n\n";
}

print "Large Contig Metrics\n";
print "Number of contigs\t", $metrics{'largeContigMetrics'}{'numberOfContigs'}, "\n";
print "Number of bases\t", $metrics{'largeContigMetrics'}{'numberOfBases'}, "\n";
print "Average contig size\t", $metrics{'largeContigMetrics'}{'avgContigSize'}, "\n";
print "N50 contig size\t", $metrics{'largeContigMetrics'}{'N50ContigSize'}, "\n";
print "Largest contig size\t", $metrics{'largeContigMetrics'}{'largestContigSize'}, "\n";
print "Q40 plus bases\t", $metrics{'largeContigMetrics'}{'Q40PlusBases'}, "\t",
	(sprintf "%.2f",	(100*$metrics{'largeContigMetrics'}{'Q40PlusBases'}/
						$metrics{'largeContigMetrics'}{'numberOfBases'})),"%\n";
print "\n";
print "All Contig Metrics\n";
print "Number of contigs\t", $metrics{'allContigMetrics'}{'numberOfContigs'}, "\n";
print "Number of bases\t", $metrics{'allContigMetrics'}{'numberOfBases'}, "\n";
print "Average contig size\t",
	(sprintf "%.0f",	$metrics{'allContigMetrics'}{'numberOfBases'}/
						$metrics{'allContigMetrics'}{'numberOfContigs'})."\n";
print "\n";

if (exists $metrics{'scaffoldMetrics'}{'numberOfScaffolds'}){
	print "Library\tPair distance average (bp)\n";
	foreach my $lib_name (sort @lib_names){
		print "$lib_name\t$metrics{$lib_name}{'pairDistanceAvg'}\n";
	}
}

