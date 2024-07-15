#!/usr/bin/env perl

use strict;

my (@Options, $verbose, $pairfile, $MPT, $thresh, $matethresh, $unknownfile, $matepairfile, $auto);

setOptions();

my %left;
my %right;
my %keys;
my %PEs;
my %MPs;
my %Unknowns;
my @dists;

my $readfile = $ARGV[1];
my $hitsfile = $ARGV[0];

unless($readfile && -r $readfile){ print "No readfile specified or not readable.\n"; &usage(); }
unless($readfile && -r $readfile){ print "No hitsfile specified or not readable.\n"; &usage(); }

if($verbose){
    print STDERR "MP_splitter.pl\n";
    print STDERR "Read file: $readfile\n";
    print STDERR "Hits file: $hitsfile\n";
    print STDERR "PE Threshold: $thresh\n";
    print STDERR "MP Threshold used: $MPT\n";
    print STDERR "MP Threshold: $matethresh\n";
    print STDERR "Verbosity: $verbose\n";
}

print STDERR "Loading hits file $hitsfile\n" if $verbose;
loadShrimpHits($hitsfile);

#look for PE and MP in the list of keys..

my $singletonCount = 0;
my $mpCount = 0;
my $peCount = 0;
my $pairCount = 0;
my $unknownpairCount = 0;

print STDERR "Looking for pairs and categorising them\n" if $verbose;

foreach my $key(keys %keys){
	#if there are less than 2 or more than 2 hits for a read pair then discard and put this read set into the unknowns..
	if($keys{$key} != 2){
		$singletonCount += $keys{$key} if $keys{$key} <= 2;
		$singletonCount += 2 if $keys{$key} > 2;
		$Unknowns{$key} ++;
		next;
	}
	if($left{$key}->{'start'} && $right{$key}->{'start'}){
    		$pairCount ++;
		#check if they are on the same ref...
		if($left{$key}->{'ref'} ne $right{$key}->{'ref'}){
			$Unknowns{$key} ++;
			$unknownpairCount += 2;
			next;
		}
	        #check if they are in opposite directions...
		if($left{$key}->{'orient'} ne $right{$key}->{'orient'}){
			#now categorise them PE, UNKNOWN or MP?
			#and put the gap distances into @dists
				push @dists, abs($left{$key}->{'start'} - $right{$key}->{'start'});
	        	if(abs($left{$key}->{'start'} - $right{$key}->{'start'}) < $thresh){
        			$peCount ++;
				$PEs{$key} ++;
	        	}
	        	elsif($MPT && abs($left{$key}->{'start'} - $right{$key}->{'start'}) < $matethresh){
        			$unknownpairCount ++;
	        	        $Unknowns{$key} ++;
		        }
        		else {
        			$mpCount ++;
            			$MPs{$key} ++;
        		}
		}
		else {
			#one of them is facing the wrong way! put them into the unknowns
			$Unknowns{$key} ++;
			$unknownpairCount += 2;
			next;
		}
    	}
    	else {
    		$singletonCount ++;
        	$Unknowns{$key} ++;
    	}
}
if($verbose){
	print STDERR "Pairs where both ends hit reference: $pairCount\n";
	print STDERR "Pairs where either but not both ends hit reference $singletonCount\n";
	print STDERR "PE count: $peCount\n";
	print STDERR "MP count: $mpCount\n";
	print STDERR "Unknown pairs: $unknownpairCount\n";
}

print STDERR "Separating read pairs based on categories and writing output\n" if $verbose;

open OUTP, ">$pairfile" or die "Couldn't open $pairfile for writing, check permissions.\n$!";
open OUTU, ">$unknownfile" or die "Couldn't open $unknownfile for writing, check permissions.\n$!";
open OUTM, ">$matepairfile" or die "Couldn't open $matepairfile for writing, check permissions.\n$!";
open OUTD, ">distancesFile.txt";

foreach (@dists){
	print OUTD "$_\n";
}
close OUTD;

open IN, $readfile or die "Couldn't open $readfile for reading.\n$!";

while(<IN>){
	chomp;
	if(/^>/){
		my $id = $_;
		my $r = $id;
		$r =~ s/\/\d$//;
		my $seq = <IN>;
		if($PEs{$r}){
			print OUTP "$id\n$seq";
		}
		elsif($Unknowns{$r}){
			print OUTU "$id\n$seq";
		}
		else{
			print OUTM "$id\n$seq";
		}
	}
}



#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"x|UseMPThresh!",  VAR=>\$MPT, DEFAULT=>1, DESC=>"Use the mate pair minimum separation threshold to call mate paired reads."},
    {OPT=>"t|pairthreshold=i",  VAR=>\$thresh, DEFAULT=>500, DESC=>"The maximum separation distance to call a paired end read."},
    {OPT=>"d|matethreshold=i",  VAR=>\$matethresh, DEFAULT=>2000, DESC=>"The minimum separation distance to call a mate pair read."},
    {OPT=>"p|pairedfile=s", VAR=>\$pairfile, DEFAULT=>"/dev/null", DESC=>"The filename to put the paired end reads into."},
    {OPT=>"u|unknownfile=s", VAR=>\$unknownfile, DEFAULT=>"/dev/null", DESC=>"The file that all the uncategorised reads will be put into."},
    {OPT=>"m|matepairfile=s", VAR=>\$matepairfile, DEFAULT=>"/dev/null", DESC=>"The file that all mate pair reads will be put into."},
    {OPT=>"a|auto!", VAR=>\$auto, DEFAULT=>0, DESC=>"Automatically creates sensibly named output files. MP_reads.fa, PE_reads.fa & Unknown_reads.fa"},
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
  
  if($auto){
	  $matepairfile = "MP_reads.fa";
	  $pairfile = "PE_reads.fa";
	  $unknownfile = "Unknown_reads.fa";
  }
  
}

sub usage {
  print "Usage: $0 [options] shrimp_hits.txt interleaved_reads.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
#	Other routines

sub loadShrimpHits {
	my $hitsfile = shift;
    
    open IN, $hitsfile or die "Couldn't open $hitsfile for reading!\n$!";
    
    while(<IN>){
    	chomp;
        next unless(/^>/);
        my @tmp = split /\t/, $_;
        my $read = $tmp[0];
		$read =~ s/\/\d$//;
		$keys{$read} ++;
		if($tmp[0] =~ /\/1$/){
			$left{$read}->{'start'} = $tmp[3];
			$left{$read}->{'ref'} = $tmp[1];
			$left{$read}->{'orient'} = $tmp[2];
	
		}
		if($tmp[0] =~ /\/2$/){
			$right{$read}->{'start'} = $tmp[3];
			$right{$read}->{'ref'} = $tmp[1];
			$right{$read}->{'orient'} = $tmp[2];
		}
   	}
    close IN;
}

#----------------------------------------------------------------------
