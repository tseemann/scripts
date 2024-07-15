#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
use Bio::Seq;

my $reffile = $ARGV[0];
my $contfile = $ARGV[1];
my $coordsfile = $ARGV[2];

unless($reffile && $contfile && $coordsfile){ die "$0 <reference fasta> <final contigs fasta> <nucmer coords>\n"; }

print "Starting fixing procedures..\nLoading coords file\n";

unless(-d "temp_data"){`mkdir temp_data`;}

my %coords;
my %scaffs;
my %contigs;
my %newscaffs;
my $blatpath = '/bio/sw/blat/';
my $revcount = 0;
my $fwdcount = 0;

open OUT, ">changelog.txt";

my $numcoords = 0;

sub dupScaffs{
	foreach my $key(keys %scaffs){
		$newscaffs{$key} = $scaffs{$key}->seq();
	}
	return 1;
}

sub resolveBothSame {
	my $l = shift;
	my $r = shift;
	
	my $lenscaff = $coords{$r}->{'send'} - $coords{$l}->{'sstart'} + 1;
	print "Length of match on scaff = $lenscaff\n";
	my $lencon = $coords{$l}->{'dist'} + $coords{$l}->{'mlen'} + $coords{$r}->{'mlen'};
	print "Length of match on contig = $lencon\n";
	
	
	my $match = $coords{$l}->{'contig'} . '$';
	unless(-e "temp_data/" .$coords{$l}->{'contig'} . ".fna"){ `fa-extract-few.pl $match < $contfile > temp_data/$coords{$l}->{'contig'}.fna`;}
	unless(-e "temp_data/" . $coords{$l}->{'scaff'} . ".fna"){ `fa-extract-few.pl $coords{$l}->{'scaff'} < $reffile > temp_data/$coords{$l}->{'scaff'}.fna`;}
	my $foo = "temp_data/$coords{$l}->{'scaff'}" . "-" . "$coords{$l}->{'contig'}";
	#`nucmer -p $foo temp_data/$coords{$l}->{'scaff'}.fna temp_data/$coords{$l}->{'contig'}.fna`;
	my $score = int(0.85 * (150 * 2 + $coords{$l}->{'dist'}));
	print "Running blat on $foo\n";
	my $x = `blat -minScore=$score -noHead -extendThroughN temp_data/$coords{$l}->{'scaff'}.fna temp_data/$coords{$l}->{'contig'}.fna -ooc=$blatpath/11.ooc /dev/stdout`;
	
	#check direction!
	
	my @t = split /\t+/, $x;
	my $cs = $t[11]+1;
	my $ce = $t[12];
	my $ss = $t[15]+1;
	my $se = $t[16];
	my $sl = $t[14];
	my ($upper, $lower, $dirn);
	
	if ($ss > $se){
		$upper = $ss;
		$lower = $se;
		$dirn = "R";
		$revcount ++;
	} else {
		$upper = $se;
		$lower = $ss;
		$dirn = "F";
		$fwdcount ++;
	}
	
	#check to see if there are any Ns in the scaff..
	
	
	my $subscaff = $scaffs{$coords{$l}->{'scaff'}}->subseq($lower, $upper);
	
	if($subscaff =~ /N/i && (($ce - $cs) <= ($se - $ss))){
		print "Scaffold contains Ns in this region. Using contig/blat to fix\n";
		#do some sanity checks.  Count the number of N's in the subscaff.. Check that the contig won't delete too many bases of non-Ns.
		#(won't delete more than 10% of the bases...)
		my $numNs = ($subscaff =~ tr/Nn/Nn/);
		my $sanity = 0;
		my $scaffLenNoN = length($subscaff) - $numNs;
		my $conrep = $contigs{$coords{$l}->{'contig'}}->subseq($cs, $ce);
		my $conNumNs = ($conrep =~ tr/Nn/Nn/);
		my $conLenNoN = length($conrep) - $conNumNs;
		unless($conLenNoN < $scaffLenNoN * 0.9){
			$sanity = 1;
		}
		if($dirn eq 'F' && $sanity){
			my $newsc;
			#got the forwards matches!
			#look at the blat output
			print $x . "\n";
			print OUT "Changing " . $coords{$l}->{'scaff'} . " using " . $coords{$l}->{'contig'} . "\n";		
			
			print "cs: $cs\tce: $ce\tss: $ss\tse: $se\tsl: $sl\tasl: " . length($scaffs{$coords{$l}->{'scaff'}}->seq()) . "\tLen conrep: " . length($conrep) . "\n";
			my $sstart = "";
			if($ss > 0){
				$sstart = $scaffs{$coords{$l}->{'scaff'}}->subseq(1 , $ss);
				print "Length sstart: " . length($sstart) . "\n";
			}
			my $send = "";
			if($se < $sl){
				$send = $scaffs{$coords{$l}->{'scaff'}}->subseq($se + 1 ,$sl);
			}
			$newsc = $sstart . $conrep . $send;
			my $newscaffobj = Bio::Seq->new(-id => $coords{$l}->{'scaff'}, -seq => $newsc);
			$scaffs{$coords{$l}->{'scaff'}} = $newscaffobj;
			my $seqout = Bio::SeqIO->new(-file => (">temp_data/" . $coords{$l}->{'scaff'} . ".fna"), -format => 'Fasta');
			$seqout->write_seq($newscaffobj);
			sleep(1);
		}
		else {
			print "Failed sanity checks. " . $contigs{$l}->{'contig'} . "\n";
		}
		
	}
	else {
		print "Scaffold does not contain Ns in this region or contig length is longer than the scaffold length.\n";
	}
	
	
	return 1;
}

print "Loading contigs\n";

#load the contigs
my $cio = Bio::SeqIO->new(-file => $contfile, -format => 'Fasta');
while(my $c = $cio->next_seq()){
	$contigs{$c->id()} = $c;
}

print "Loading reference\n";
#load the scaffs
my $sio = Bio::SeqIO->new(-file => $reffile, -format => 'Fasta');
while(my $s = $sio->next_seq()){
	$scaffs{$s->id()} = $s;
}

&dupScaffs();


open IN, "$coordsfile";

while(<IN>){
	chomp;
	my @t = split;
	my @n = split /_/, $t[12];
	my $name = $n[0] . "_" . $n[1] . "_" . $n[2];
	my $contig = $n[0] . "_" . $n[1];
	unless($t[3] < $t[2]){
		$coords{$name}->{'contig'} = $contig;
		$coords{$name}->{'dir'} = $n[2];
		$coords{$name}->{'dist'} = $n[3];
		$coords{$name}->{'scaff'} = $t[11];
		$coords{$name}->{'slength'} = $t[7];
		$coords{$name}->{'sstart'} = $t[0];
		$coords{$name}->{'send'} = $t[1];
		$coords{$name}->{'cstart'} = $t[2];
		$coords{$name}->{'cend'} = $t[3];
		$coords{$name}->{'idty'} = $t[6];
		$coords{$name}->{'mlen'} = $t[5];
		$numcoords ++;
	}
}

print "$numcoords nucmer coords loaded.\nWorking through contigs\n";

foreach my $key (sort keys %coords){
	
	if($key =~ /left/){
		print "Contig: " . $coords{$key}->{'contig'} . "\n";
		my $right = $coords{$key}->{'contig'} . "_right";
		if($coords{$right}){
			if($coords{$key}->{'scaff'} eq $coords{$right}->{'scaff'}){
				print "Both ends are on the same scaffold!- " . $coords{$key}->{'scaff'} . "\n";
				&resolveBothSame($key, $right);			
			}
			else {
				print "Ends are on different scaffolds!- left: " . $coords{$key}->{'scaff'} . " right: " . $coords{$right}->{'scaff'} . "\n";
			}
		}
		else {
			print "Left: " . $coords{$key}->{'scaff'} . " Right hand end is not on any scaffold!\n";
		}
	}
	elsif(($key =~ /right/) && !($coords{$coords{$key}->{'contig'} . "_left"})){
		print "Contig: " . $coords{$key}->{'contig'} . "\n";
		print "Right: " . $coords{$key}->{'scaff'} . " Left hand end is not on any scaffold!\n";
	}
}

print "Number of forward contig matches: $fwdcount\n";
print "Number of reverse contig matches: $revcount\n";
