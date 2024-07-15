#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastUtils qw(read_fasta read_fasta_Bio read_fastq read_fastq_Bio write_fastq write_fasta assert detect_filetype);
use Bio::SeqIO;

my(@Options, $verbose, $interleave, $fasta, $newid, $nons, $trunc, $nopoly, $bio);
setOptions();

assert(-r $ARGV[0], "can not open left reads: $ARGV[0]");
if ($interleave) {
	assert(@ARGV >= 2, "please provide two .fq files as input");
	assert(-r $ARGV[1], "can not open right reads: $ARGV[1]");
}

#
#
#	Find out what type of file(s) we have been given and decide what to
#	do with them...
#
#

my $ltype = &detect_filetype($ARGV[0]);
my $rtype = &detect_filetype($ARGV[1]) if $interleave;

#print STDERR "LType: $ltype\n";
#print STDERR "RType: $rtype\n" if $interleave;

if($ltype eq "Unknown") { die "Unknown file type specified: $ARGV[0]\n"; }
if ($interleave){
	if($rtype eq "Unknown") { die "Unknown file type specified: $ARGV[1]\n"; }
}

my $fastain = 0;

if($ltype =~ /Fasta/){
	$fastain = 1;
	#Since with fasta in we can't specify fastq out, make sure that we don't try...
	$fasta = 1;
}
if($ltype =~ /Bio$/){
	$bio = 1;
}
#now if there are two files, they need to be the same format!
if($interleave){
	unless($ltype eq $rtype) {
		die "The files $ARGV[0] and $ARGV[1] must have the same file format, they are $ltype and $rtype respectively..\n";
	}
	if($rtype =~ /Bio$/){
		$bio = 1;	#even if one is bio easier to make both bio...
	}
}

#
#
# now choose the appropriate reader and writer functions
#
#

my $writer_func = $fasta ? \&write_fasta : \&write_fastq;
my $reader_func;
if($bio){
	$reader_func = $fastain ? \&read_fasta_Bio : \&read_fastq_Bio;
}
else {
	$reader_func = $fastain ? \&read_fasta : \&read_fastq;
}

#
#
#	open the file(s) in different manners depending on $bio...
#
#

my ($L, $R);

#print STDERR "Bio = $bio\tfasta_out = $fasta\t fasta in = $fastain\tARGV[0] = $ARGV[0]\n";

if($bio){
	$L = Bio::SeqIO->newFh(-file => $ARGV[0], -format => $ltype);
	$R = Bio::SeqIO->newFh(-file => $ARGV[1], -format => $rtype) if $interleave;
}
else {
	open $L, '<', $ARGV[0];
	open $R, '<', $ARGV[1] if $interleave;
}

#
#
#	Write some info to STDERR for users who ask for it!
#
#

if ($verbose) {
	print STDERR "Interleaving $ARGV[0] + $ARGV[1]\n" if $interleave;
	print STDERR "Using Bio::Perl to read the files\n" if $bio;
	print STDERR "Reading fast", ($fastain ? 'a' : 'q'), ".\n";
	print STDERR "Removing pairs with 'N's.\n" if $nons;
	print STDERR "Removing monopolymers.\n" if $nopoly;
	print STDERR "Truncating reads to $trunc bp.\n" if $trunc;
	print STDERR "Output is FAST", ($fasta ? 'A' : 'Q'), ".\n";
	print STDERR "Please be patient...\n";
}

#
#
#	Start doing the business!
#
#

my $nread=0;
my $nwrote=0;

while ($bio ? 1 : not eof $L) {
	
	my $l = $reader_func->($L);
	unless($l) { last; }
	my $r = $reader_func->($R) if $interleave;
	$nread++;

	# truncate first
	if ($trunc > 0) {
		$l->[1] = substr($l->[1], 0, $trunc);
		$l->[2] = substr($l->[2], 0, $trunc) unless $fastain;
		if ($interleave){
			$r->[1] = substr($r->[1], 0, $trunc);
			$r->[2] = substr($r->[2], 0, $trunc) unless $fastain;
		}
	}

	# then skip read if either left or right has N in it
	
	if($interleave){
		next if $nons and index($l->[1].$r->[1], 'N') >= $[ ;
	}
	else {
		next if $nons and index($l->[1], 'N') >= $[ ;
	}
  
	# don't allow monopolymers through ?
	if($interleave){
		next if $nopoly and ( $l->[1] =~ m/^(.)\1*$/i or $r->[1] =~ m/^(.)\1*$/i );
	}
	else {
		next if $nopoly and ( $l->[1] =~ m/^(.)\1*$/i);
	}
                       
	$nwrote++;
	
	#change the ids if required
	if ($newid) {
		$l->[0] = "$nwrote/1";
		$r->[0] = "$nwrote/2" if $interleave;
	}

	$writer_func->(\*STDOUT, $l);
	$writer_func->(\*STDOUT, $r) if $interleave;
}
	
printf STDERR "Parsed %d ", $nread;
$interleave ? print STDERR "Pairs." : print STDERR "Reads.";

printf STDERR " Discarded %d. Wrote %d (%.2f%%)\n", 
	$nread-$nwrote, 
	$nwrote, ($nwrote*100/$nread);
	
print STDERR "Script took " . (time() - $^T) . " seconds to run.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"i|interleave!", VAR=>\$interleave, DEFAULT=>0, DESC=>"Interleave two paired end files"}, 
		{OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
		{OPT=>"n|nons!", VAR=>\$nons, DEFAULT=>0, DESC=>"Skip pairs with Ns"},
		{OPT=>"p|nopoly!", VAR=>\$nopoly, DEFAULT=>0, DESC=>"Skip monopolymer reads"},
		{OPT=>"t|trunc=i", VAR=>\$trunc, DEFAULT=>0, DESC=>"Truncate reads to this length (0 = no truncation)"},
		{OPT=>"d|newid!", VAR=>\$newid, DEFAULT=>0, DESC=>"Generate new, compact IDs"},
		{OPT=>"b|bio!", VAR=>\$bio, DEFAULT=>0, DESC=>"Force the use of the Bio::SeqIO module to read sequences"},
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
	print "Usage: $0 [options] left.fq [right.fq] > interleaved.fa\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
