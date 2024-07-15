#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";

my(@Options, $verbose, $lp, $rp, $maxpairsep);
setOptions();

my $file = $ARGV[0];

unless(open IN, $file){ die "Unable to open shrimp file $file for reading.\n$!"; }

my %reads;
my %LHS;
my %RHS;
my $singles = 0;
my $pairs = 0;
my $halfpairs = 0;
my $count = 0;

while(<IN>){
	unless(/^>/){ next; }
	chomp;
	my @tmp = split /\t/, $_;
	my ($q, $ref, $str, $cs, $ce, $rs, $re, $rl, $score, $estring, $rseq) = @tmp;
	unless($rseq) { die "Re-run nesoni with '--shrimp-options -R' and start again!\n"; }
	
	#get rid of the >
	$q =~ s/^>//;
	
	#got read, is pair?
	if($q =~ m/($lp|$rp)$/){
		#pair!
		#have we seen it?
		if($reads{$q}->{'num'}){
			$reads{$q}->{'num'}++;
			#is the score the same as the previous?
			if($score >= $reads{$q}->{'score'}){
				$q .= "." . $reads{$q}->{'num'}++;
				$singles ++;
				&processSingle($q, $ref, $str, $cs, $ce, $rs, $re, $rl, $score, $estring, $rseq, 1);
			}
			else {
				next;
			}
		}
		else {
			#haven't seen it!
			$reads{$q}->{'num'} = 1;
			$reads{$q}->{'score'} = $score;
			$singles ++;
			&processSingle($q, $ref, $str, $cs, $ce, $rs, $re, $rl, $score, $estring, $rseq);
			#my $key = $q;
			#$key = s/($lp|$rp)$//;
			#print STDERR "$key\n";
			
			#is it LHS or RHS?
			#if($q =~ m/$lp$/){
				#its left!
				#store left
			#	$LHS{$key} = [$q, $ref, $str, $cs, $ce, $rs, $re, $rl, $score, $estring, $rseq];
			#}
			#else {
				#its right!
				#store right
			#	$RHS{$key} = [$q, $ref, $str, $cs, $ce, $rs, $re, $rl, $score, $estring, $rseq];
			#}
			#do we have both?
			#if($LHS{$key} && $RHS{$key}){
				#we have both!
			#	$pairs += 2;
			#	&processPair($LHS{$key}, $RHS{$key});
			#}
		}
	}
	else {
		#no pair!
		#have we seen it?
		if($reads{$q}->{'num'}){
			$reads{$q}->{'num'}++;
			#is the score the same as the previous?
			if($score >= $reads{$q}->{'score'}){
				$q .= "." . $reads{$q}->{'num'}++;
				$singles ++;
				&processSingle($q, $ref, $str, $cs, $ce, $rs, $re, $rl, $score, $estring, $rseq, 1);
			}
			else {
				next;
			}
		}
		else {
			#haven't seen it!
			$reads{$q}->{'num'} = 1;
			$reads{$q}->{'score'} = $score;
			$singles ++;
			&processSingle($q, $ref, $str, $cs, $ce, $rs, $re, $rl, $score, $estring, $rseq);
		}
	}
	
	$count ++;
	if($count % 100000 == 0){
		print STDERR "\rProcessed $count reads";
	}
	
}

print "\n";

#need to look for things in LHS that don't exist in RHS and vice-versa..
#and process them as half pairs...  (bloody SAM and its bit flags..)
foreach my $key (keys %LHS){
	unless($RHS{$key}){
		$halfpairs++;
		&processHalfPair($LHS{$key});
	}
}

foreach my $key (keys %RHS){
	unless($LHS{$key}){
		$halfpairs++;
		&processHalfPair($RHS{$key});
	}
}


print STDERR "Singles processed: $singles\n";
print STDERR "Pairs processed: $pairs\n";
print STDERR "Half pairs processed: $halfpairs\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"l|left-pair-suffix", VAR=>\$lp, DEFAULT=>"/1", DESC=>"The left hand suffix of the read name"},
		{OPT=>"r|right-pair-suffix", VAR=>\$rp, DEFAULT=>"/2", DESC=>"The right hand suffix of the read name"},
		{OPT=>"s|maxpairsep=i", VAR=>\$maxpairsep, DEFAULT=>500, DESC=>"The maximum pair separation for the reads to be counted as paired."},
	);

	

	&GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

	(@ARGV < 1) && (usage());

	# Now setup default values.
	foreach (@Options) {
		if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
		${$_->{VAR}} = $_->{DEFAULT};
		}
	}
}

sub usage {
	print "Usage: $0 [options] shrimp_hits.txt > \n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------

#editString2CIGAR
#ripped off from Nils Homer's shrimp2sam.pl supplied with SHRiMP in the utils directory.
sub editString2CIGAR {
	
	my $es = shift;
	my $rs = shift;	#used for soft clipping of reads...
	my $re = shift;	#used for soft clipping of reads...
	my $rl = shift;	#used for soft clipping of reads...
	
	my $tmp = uc($es);
	 
	my $out = ""; #$tmp . "\t";
	
	my $prev_m = 0;
	
	while(0 < length($tmp)) {
		if($tmp =~ m/^(\d+)/) { # ref match
			$prev_m += $1;
			$tmp = substr($tmp, length($1));
		}
		elsif($tmp =~ m/^([ACGTN])/) { # snp
			$prev_m++;
			$tmp = substr($tmp, 1);
		}
		elsif($tmp =~ m/^X/) { # crossover
			$tmp = substr($tmp, 1);
		}
		else { # indel 
			if(0 < $prev_m) {
				$out .= "$prev_m"."M"; $prev_m=0;
			}
			if($tmp =~ m/^\(([ACGTN]+)\)/) { # insertion
				$out .= "".length($1)."I";
				$tmp = substr($tmp, length($1) + 2);
			}
			elsif($tmp =~ m/^(\-+)/) { # deletion
				$out .= "".length($1)."D";
				$tmp = substr($tmp, length($1));
			}
			else {
				print "$tmp\n";
				die;
			}
		}
	}
	if(0 < $prev_m) {
		$out .= "$prev_m"."M"; $prev_m=0;
	}
	unless($rs == 1){
		my $n = $rs -1;
		my $t = $out;
		$out = $n . "S" . $t;
	}
	unless($re == $rl){
		my $n = $rl - $re;
		my $t = $out;
		$out = $t . $n . "S";
	}
	
	return $out;
}

sub revcomp {
	my $s = shift;
	$s =~ tr/atcgATCG/tagcTAGC/;
	my $out = reverse $s;
	return $out;
	
}

sub processSingle {
	my ($q, $ref, $str, $cs, $ce, $rs, $re, $rl, $score, $estring, $rseq, $prim) = @_;
	my $out = "";
	
	#set the query..
	$out .= $q . "\t";
	
	#build the flag from the only possible values for a single...
	my $flag = 0x0000;
	$flag |= 0x0010 if $str eq '-';
	$flag |= 0x0100 if $prim;
	$out .= $flag . "\t";
	
	#add the reference..
	$out .= "$ref\t";
	
	#add the read contig start point..
	$out .= "$cs\t";
	
	#add the mapping quality (for Shrimp its unknown so 255...)
	$out .= "255\t";
	
	#convert the edit string to CIGAR format and add it to the line...
	$out .= &editString2CIGAR($estring, $rs, $re, $rl) . "\t";
	
	#MRNM (Mate reference sequence name..), MPOS and ISize are * 0 0 respectivley as these are not paired reads.
	$out .= "*\t0\t0\t";
	
	#read sequence..
	$out .= uc($rseq) . "\t" if $str eq '+';
	$out .= uc(revcomp($rseq)) . "\t" if $str eq '-';
	
	#read quality..  Not present in SHRiMP data...
	$out .= "*\t";
	
	#TAGs..  Only one needed here..  The alignment score tag generated by SHRiMP..
	$out .= "AS:i:$score";
	
	$out .= "\n";
	
	print $out;
	return 1;
}

sub processPair {
	
	my $lhs = shift;
	my $rhs = shift;
	
	my $lout = "";
	my $rout = "";
	
	#left hand side...
	my ($lq, $lref, $lstr, $lcs, $lce, $lrs, $lre, $lrl, $lscore, $lestring, $lrseq) = @$lhs;
	#right hand side...
	my ($rq, $rref, $rstr, $rcs, $rce, $rrs, $rre, $rrl, $rscore, $restring, $rrseq) = @$rhs;
	
	print abs($rcs - $lcs) . "\n";
	
	my $str = join "\t", @$lhs;
	my $rstr = join "\t", @$rhs;
	
	print "$str\n$rstr\n";
	
	if(abs($rcs - $lcs) > $maxpairsep){
		
		$pairs -= 2;
		$halfpairs +=2;
		&processHalfPair($lhs);
		&processHalfPair($rhs);
		return 1;
	}
	
	
	#left hand side first!
	#set the queries..
	$lout .= "$lq\t"; $rout .= "$rq\t";
	
	#now for the flags!
	#left hand side first!
	my $flag = 0x0000;
	$flag |= 0x0001; #cause its a paired read!
	
	
	return 1;
}

sub processHalfPair {
	return 1;
}
