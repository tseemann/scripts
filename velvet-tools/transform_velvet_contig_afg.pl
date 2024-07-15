#!/usr/bin/perl
#
#
#	A script to convert afg file to their associated names
#	Author: Jason Tsai jit@sanger.ac.uk
#	Last modified: 02.06.2009
#
#

if (@ARGV != 3) {
	print "transform_velvet_contig_afg.pl original_afg input_fasta contigs.fa\n\n" ;
	exit ;
}


$filenameA = $ARGV[0];
$filenameB = $ARGV[1];
$filenameC = $ARGV[2];




open $FILEA, "< $filenameA";
open $FASTA, "< $filenameB";
open $FASTA2, "< $filenameC";

# a hack to convert velvet output in afg
my $iid = 1 ;
my %correspond ;

while (<$FASTA>) {
	if (/^>(\S+)/) {
		$correspond{$iid} = "$1" ;
		$iid++ ;
	
	}
}
close(FASTA) ;


my $iid = 1 ;
my %contig_correspond ;

while (<$FASTA2>) {
	if (/^>(\S+)/) {
		$contig_correspond{$iid} = "$1" ;
		$iid++ ;
	
	}
}
close(FASTA2) ;



# a fake id for the fragment size
$iid = 1 ;

while(<$FILEA>) {
	if (/FRG/) {
		my $line = $_ ;
		print "$line" ;
		
		$line = <$FILEA> ;
		print "$line" ;
		
		$line = <$FILEA> ;
		print "$line" ;

		print "iid:$iid\n" ;
		$iid++ ;

		print "typ:I\n" ;

		$line = <$FILEA> ;
		print "$line" ;

	}
	elsif (/RED/) {
		my $line = $_ ;
		print "$line" ;
		
		$line =  <$FILEA> ;
		chomp($line) ;
		@lines = split /:/, $line ;
		print "$line\n" ;
		
		$line =  <$FILEA> ;
		print "eid:" . $correspond{$lines[1]} . "\n" ;		

	}
	elsif (/\{CTG/) {
		my $line = $_ ;
		print "$line" ;
		
		$line =  <$FILEA> ;
		chomp($line) ;
		@lines = split /:/, $line ;
		print "$line\n" ;
		
		$line =  <$FILEA> ;
		print "eid:" . $contig_correspond{$lines[1]} . "\n" ;		

	}
	else {
		print $_ ;
	}

	

	
}
