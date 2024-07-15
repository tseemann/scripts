#!/usr/bin/perl -w

use strict;

my $reffile = $ARGV[0];
my $confile = $ARGV[1];
my $endsize = 300;
my $mincontig = 1000;

my $usage = "$0 <reference file> <contigs file>\n";

unless($reffile && -e $reffile && -r $reffile){die "Problem with ref file!\n$usage\n"; }
unless($confile && -e $confile && -r $confile){die "Problem with contig file!\n$usage\n"; }

#first make the names of the contigs more user friendly...
print "Cleaning contig names..\n";
`fa-clean_velvet_contig_names.pl $confile clean_contigs.fna`;
sleep(5);

print "Beginning pre-processing..\nGetting large contigs..\n";
`fa-filter-size.pl -min $mincontig clean_contigs.fna > clean_contigs_500.fna`;
sleep(5);

print "Re-orienting the contigs to match reference..\n";
`fa-reorient.pl -r $reffile -i clean_contigs_500.fna > clean_contigs_500.fna.ro`;
sleep(5);

print "Producing fasta of contig ends..\n";
`fa-extract_contig_ends.pl -t 0 -b $endsize clean_contigs_500.fna.ro > Ends_clean_contigs_500.fna.ro`;
sleep(5);

print "'Nucmer'ing the ends against the reference..\n";
`nucmer $reffile Ends_clean_contigs_500.fna.ro`;
sleep(5);

print "Producing coordinates file..\n";
`show-coords -clTH out.delta > nucmer-coords.txt`;

print "Making output directory and splitting scaffolds file..\n";
`mkdir temp_data`;
chdir "temp_data";
`fa-split.pl --template=%i.fna ../$reffile`;

print "Finished pre-processing and end matching!\n";
