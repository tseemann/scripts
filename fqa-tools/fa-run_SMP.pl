#!/usr/bin/env perl
#       fa-run_SMP.pl
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
use File::Temp qw (tempfile tempdir);
use threads;
use threads::shared;
use Carp;

my(@Options, $verbose, $command, $max_threads, $num_seqs, $dryrun);
setOptions();

my $num_threads : shared = 0;

my $fastafile = $ARGV[0];

unless($fastafile && -r $fastafile) { print "No fasta file supplied or not readable..\n"; usage(); }

#make a temporary directory to blat files into....

my $dir = tempdir( CLEANUP => 1);

#make the temp split up fasta files...

my $in = Bio::SeqIO->new(-file => $fastafile, -format => 'Fasta');
my $fcount = 0;
my $out;
my $scount = $num_seqs +1;
my $outname;
my @fastafilenames;
my @outfilename;

while (my $seq = $in->next_seq) {
  if ($scount >= $num_seqs) {
    $fcount++;
    $scount=0;
    $outname = "$fcount" . ".fa";
    print STDERR "Splitting to: $outname\n" if $verbose;
    $out = Bio::SeqIO->new(-file => ">$dir/$outname", '-format' => 'Fasta')
      or die "coult not write to '$dir/$outname'";
    push @fastafilenames, $outname;
  }
  $out->write_seq($seq);
  $scount++;
}

print STDERR "Split $fastafile into " . scalar @fastafilenames . " files\n" if $verbose;

my @threads;

#ok so now we want to make a loop to loop through all the files created and pass them one by one to the threaded runCommand routine..
foreach my $fn(@fastafilenames){
	my $outfile = $fn;
	$outfile =~ s/.fa//;
	my $thr_num = $outfile;
	$outfile .= ".out";
	push @outfilename, $outfile;
	my $stdoutfile = $outfile . ".stdout";
	my $stderrfile = $outfile . ".stderr";
	
	my $comm = $command;
	$comm =~ s/%i/$dir\/$fn/g;
	if($comm =~ m/>(\s+)?%o/){
		$comm =~ s/%o/$dir\/$outfile/g;
	}
	elsif($comm =~ m/%o/){
		$comm =~ s/%o/$dir\/$outfile/g;
		$comm .= " 1> $dir/$stdoutfile";
	}
	else {
		$comm .= " > $dir/$outfile";
	}
	
	$comm .= " 2> $dir/$stderrfile";
	
	#print "Starting to run $comm\n";
	while($num_threads >= $max_threads){
		sleep (1);
	}
	$threads[$thr_num] = threads->create(\&runCommand, $comm, $thr_num);
	{
		lock($num_threads);
		$num_threads ++;
	}
	
}

for my $thr (threads->list) {
	$thr->join;
}

#now go through the output files and cat them and spit them out...
unless($dryrun){
	foreach my $fn (@outfilename){
		#print $fn . "\n";
		system("cat $dir/$fn") if(-r "$dir/$fn");
	}
}

unless($dryrun){
	foreach my $fn (@outfilename){
		#print $fn . "\n";
		system("cat $dir/$fn.stdout >> $fastafile.out.stdout") if(-r "$dir/$fn.stdout");
	}
}

unless($dryrun){
	foreach my $fn (@outfilename){
		#print $fn . "\n";
		system("cat $dir/$fn.stderr >> $fastafile.out.stderr") if(-r "$dir/$fn.stderr");
	}
}


#----------------------------------------------------------------------
# Run the command subroutine - (its the thread rountine!)

sub runCommand {
	my $comm = shift;
	my $thread_name = shift;
	
	print STDERR "Running $comm\n" if $verbose;
	if($dryrun){
		sleep (int(rand(5)));
	}
	else{
		system($comm) == 0 or carp "system $comm failed\n$?\n";
	}
	print STDERR "Finished running $comm\n" if $verbose; 
	{
		lock($num_threads);
		$num_threads --;
	}
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	my $thmax = num_cpu();

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"c|command=s", VAR=>\$command, DEFAULT=>"", DESC=>"The command line to run*"},
		{OPT=>"t|threads=i", VAR=>\$max_threads, DEFAULT=>$thmax, DESC=>"The number of threads to run at once."},
		{OPT=>"n|numseqs=i", VAR=>\$num_seqs, DEFAULT=>100, DESC=>"The number of sequences to put into each split file."},
		{OPT=>"d|dryrun!", VAR=>\$dryrun, DEFAULT=>0, DESC=>"Perform a dry run without running actual commands."},
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
	print "Usage: $0 [options] fastafile.fa > outputfile\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	
	print "\nFor the -c|command option, use a normal fastafile utility command in \"\"s.\nUse:\n\t%i as placemark for the input file and\n\t%o as placemark for the output file.\n";
	print "Example:\n\tTo run blastp on a fasta file in SMP mode using 8 threads/cpus, use the following:\n";
	print "\t$0 -c \"blastp -db database -query %i -evalue 0.1 -out %o\" -t 8 my_fasta_file_of_proteins.faa\n";
	
	exit(1);
}
 
#----------------------------------------------------------------------


# 	num_cpu
#	It returns the number of cpus present in the system if linux.
#	If it is MAC then it returns the number of cores present.
#	If the OS is not linux or Mac then it returns 1.
#	Written by Torsten Seemann 2009 (linux) and Mikael Brandstrom Durling 2009 (Mac).

sub num_cpu {
    if ( $^O =~ m/linux/i ) {
        my ($num) = qx(grep -c ^processor /proc/cpuinfo);
        chomp $num;
        return $num if $num =~ m/^\d+/;
    }
	elsif( $^O =~ m/darwin/i){
		my ($num) = qx(system_profiler SPHardwareDataType | grep Cores);
		$num =~ /.*Cores: (\d+)/;
		$num =$1;
		return $num;
	}
    return 1;
}
