#!/usr/bin/env perl
#       phantom2.pl
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
use Bio::SimpleAlign;
use Bio::AlignIO;

my(@Options, $verbose, $scaffolds_file, $reads_file, $num_iters, $dir, $flank, $cpus, $k_val);
setOptions();

#other future options.
my $min_match = 95;	#the minimum %id for a blast match
my $min_match_length = $flank - 50;	#the minimum length for a blast match

#get the file names
$scaffolds_file = $ARGV[0];
$reads_file = $ARGV[1];

unless(-r $scaffolds_file){ die "Unable to read $scaffolds_file or it doesn't exist.\n $!\n"; }
unless(-r $reads_file){ die "Unable to read $reads_file or it doesn't exist.\n $!\n"; }

#unless(&checkRequirements()){ die "Some of the requirments could not be met. $! \n"; }

print "phantom2.pl - Copyright Simon Gladman 2010, CSIRO.\nVersion 0.1\n";

#make the output dir..
unless(-d $dir) { unless(mkdir $dir){ die "Unable to create output dir. $!\n"; } }
chdir $dir;

#copy the original scaffolds file to the out_dir...
system("cp ../$scaffolds_file .");

#start iterating...
my $iter_count = 1;

print  "\n\nPerforming $num_iters iterations!\n\n";

while($num_iters > 0){
	
	print  "Beginning the $iter_count iteration!\n\n";
	if($iter_count > 1){
		my $oic = $iter_count - 1;
		$scaffolds_file = "new_scaffolds_${oic}.fa";
		print  "New Scaffolds file = $scaffolds_file\n";
	}
	
	#load the scaffolds file into memory...
	my %scaffolds;
	my @scaff_order;
	my %joins;
	my $scaffio = Bio::SeqIO->new(-file => "$scaffolds_file", -format => 'Fasta');
	while(my $seq = $scaffio->next_seq){
		$scaffolds{$seq->id} = $seq->seq;
		push @scaff_order, $seq->id;
	}
	
	#split the scaffolds file to make contig fa's..
	my $cmd = "fa-split.pl $scaffolds_file";
	print  "$cmd\n";
	system($cmd);
	
	#get the contigs and make flanking sequences...
	$cmd = "fa-extract_contig_ends.pl -b 300 -t 0 -f $scaffolds_file > contig_flanks_${iter_count}.fa";
	print  "$cmd\n";
	system($cmd);
	
	#split the contig flanks up...
	$cmd = "fa-split.pl contig_flanks_${iter_count}.fa";
	print  $cmd . "\n";
	system($cmd);

	#next thing to do is to map the reads to these ends...
	$cmd = "/bio/sw/shrimp2/bin/gmapper -N $cpus -o 1 -h 80% ../$reads_file contig_flanks_${iter_count}.fa > shrimp_hits_${iter_count}.txt ";
	print  $cmd . "\n";
	system($cmd);
	
	#iterate through all flanking sequences and make new extensions...
	my $flank_io = Bio::SeqIO->new(-file => "contig_flanks_${iter_count}.fa", -format => "Fasta");
	
	#skip the first flanking sequence as it is the beginning and really only matches the end of the genome or plasmid for circulars..
	my $seq = $flank_io->next_seq();
	
	#then we need to work on them in pairs...
	
	my $count = 0;
	
	while(1){
		
		my ($seqL, $seqR);
		$seqL = $flank_io->next_seq() or last;
		$seqR = $flank_io->next_seq() or last;
		
		#open the log file for this gap.
		my $logfile;
		print "\n*******\n" . $seqL->id . "\n";
		$seqL->id =~ m/(\d+)_((left)|(right))_-?\d+$/;
		$logfile = "gap_" . $1 . ".txt";
		print "Logfile name = $logfile\n\n";
		if(-e $logfile){
			open LOG, ">>$logfile" or die "Couldn't open logfile for writing. $logfile\n$!\n";
		}
		else {
			open LOG, ">$logfile" or die "Couldn't open logfile for writing. $logfile\n$!\n";
		}
		my $alIO = Bio::AlignIO->new(-fh => \*LOG, -format => "clustalw");
		print LOG "Iteration $iter_count\n";
			
		print  "\n\nProcessing " . $seqL->id . " and " . $seqR->id . "\n";
		#grep out the read names...
		$cmd = "grep -F \"" . $seqL->id . "\" shrimp_hits_${iter_count}.txt | cut -f 1 | " . 'sed "s/\/[12]$//" ' . "> read_namesL.txt";
		print  "$cmd\n";
		system($cmd);
		$cmd = "grep -F \"" . $seqR->id . "\" shrimp_hits_${iter_count}.txt | cut -f 1 | " . 'sed "s/\/[12]$//" ' . "> read_namesR.txt";
		print  "$cmd\n";
		system($cmd);
		
		#grep out the reads...
		$cmd = "grep -F -f read_namesL.txt -A 1 ../$reads_file | grep -v '^--\$' > reads_for_assemblyL.fa";
		print  "$cmd\n";
		system($cmd);
		$cmd = "grep -F -f read_namesR.txt -A 1 ../$reads_file | grep -v '^--\$' > reads_for_assemblyR.fa";
		print  "$cmd\n";
		system($cmd);
		
		#cat the reads for a joint assembly.
		$cmd = "cat reads_for_assemblyL.fa reads_for_assemblyR.fa > reads_for_assembly_all.fa";
		print  "$cmd\n";
		system($cmd);
		
		#check the file size
		unless(-s "reads_for_assembly_all.fa" > 0){
			print "No reads to assemble!\n";
			print LOG "No reads to process in iteration #$iter_count\n";
			close LOG;
			next;
		}
		
		#assemble the reads
		$cmd = "velveth tempAll $k_val -shortPaired reads_for_assembly_all.fa 1> /dev/null 2> /dev/null";
		print  "$cmd\n";
		system($cmd);
		
		$cmd = "velvetg tempAll -exp_cov auto -cov_cutoff auto -min_contig_lgth 150 -scaffolding no 1> /dev/null 2> /dev/null";
		print  "$cmd\n";
		system($cmd);
		
		#move and rename the contigs
		my $x = $seqL->id;
		my $y = $seqR->id;
		$x =~ s/[\/:]/_/g;
		$y =~ s/[\/:]/_/g;
		
		$cmd = "mv tempAll/contigs.fa contigsAll_${y}.fa";
		print  "$cmd\n";
		system($cmd);
		system("rm -r tempAll");
		
		#now fa the assembly output..
		$cmd = "fa contigsAll_${y}.fa 1>&2";
		print  "$cmd\n";
		system($cmd);
		
		#cat the assembly contigs and the left/right flank seqs
		$cmd = "cat contigsAll_${y}.fa ${x}.seq ${y}.seq > forCAPall.fa";
		print  "$cmd\n";
		system($cmd);
		
		#give it all to CAP3 for secondary assembly
		$cmd = "cap3 forCAPall.fa > capOUTall.txt";
		print  "$cmd\n";
		system($cmd);
		
		#now fa the cap output..
		$cmd = "fa forCAPall.fa.cap.contigs 1>&2";
		print  "$cmd\n";
		system($cmd);
		
		#load the cap assembly contigs into a hash called queries!
		my %queries;
		my $qcount = 0;
		my $queryio = Bio::SeqIO->new(-file => "forCAPall.fa.cap.contigs", -format => 'Fasta');
		while(my $seq = $queryio->next_seq()){
			$queries{$seq->id} = $seq->seq;
			$qcount ++;
		}
		
		#get the correct contigs together and formatdb them from the split scaffolds files...
		#process file names..
		my $x_full = $x;
		$x_full =~ s/_((left)|(right))_.*$/\.seq/;
		print "x_full: $x_full\n";
		my $y_full = $y;
		$y_full =~ s/_((left)|(right))_.*$/\.seq/;
		print "y_full: $y_full\n";
		
		#formatdb the contig for the left hand side of the gap
		$cmd = "formatdb -i $x_full -p F -t left_hand -n left_hand";
		print  "$cmd\n";
		system($cmd);

		#formatdb the contig for the left hand side of the gap
		$cmd = "formatdb -i $y_full -p F -t right_hand -n right_hand";
		print  "$cmd\n";
		system($cmd);
		
		#blast the cap contigs versus the left hand end! Store the best result with a bunch of stats...
		#for display and checking purposes only I've done a blast to STDOUT here...
		$cmd = "blastn -db left_hand -query forCAPall.fa.cap.contigs -evalue 0.00001 -outfmt 6";
		print  "\n\n$cmd\n";
		system($cmd);
		
		my %lh_bl_result;
		$lh_bl_result{score} = 0;
		my $lh_sio = Bio::SearchIO->new(-file=>"blastn -db left_hand -query forCAPall.fa.cap.contigs -evalue 0.00001 |", -format=>'blast');
		while (my $res = $lh_sio->next_result) {
			next if $res->no_hits_found;
			while (my $hit = $res->next_hit) {
				while (my $hsp = $hit->next_hsp) {
					if($hsp->bits > $lh_bl_result{score}){
						#we have a better hsp!  Store it!
						$lh_bl_result{score} = $hsp->bits;
						$lh_bl_result{queryname} = $res->query_name;
						$lh_bl_result{hitname} = $hit->name;
						$lh_bl_result{id} = $hsp->percent_identity;
						$lh_bl_result{hsplength} = $hsp->length;
						$lh_bl_result{querylength} = $res->query_length;
						$lh_bl_result{hitlength} = $hit->length;
						$lh_bl_result{qstrand} = $hsp->strand("query");
						$lh_bl_result{hstrand} = $hsp->strand("hit");
						$lh_bl_result{querystart} = $hsp->query->start;
						$lh_bl_result{queryend} = $hsp->query->end;
						$lh_bl_result{subjectstart} = $hsp->sbjct->start;
						$lh_bl_result{subjectend} = $hsp->sbjct->end;
						$lh_bl_result{aln} = $hsp->get_aln;
					}
				}
			}
		}
		if($lh_bl_result{score} == 0){
			print "\nNo blast hits found on left hand end!\n";
			print LOG "\nNo blast hits found on left hand end!\n";
		} else {
			print "\nLeft hand best hit/hsp etc:\n";
			print "\tQuery\t",$lh_bl_result{queryname}, "\n";  
			print "\tHit\t", $lh_bl_result{hitname}, "\n";
			print "\tScore\t", $lh_bl_result{score}, "\n";
			print "\t%id\t", $lh_bl_result{id}, "\n";
			print "\tMLength\t", $lh_bl_result{hsplength}, "\n";
			print "\tQLength\t", $lh_bl_result{querylength}, "\n";
			print "\tSLength\t", $lh_bl_result{hitlength}, "\n";
			print "\tQStrand\t", $lh_bl_result{qstrand}, "\n";
			print "\tHStrand\t", $lh_bl_result{hstrand}, "\n";
			print "\tQStart\t", $lh_bl_result{querystart}, "\n";
			print "\tQEnd\t", $lh_bl_result{queryend}, "\n";
			print "\tSStart\t", $lh_bl_result{subjectstart}, "\n";
			print "\tSEnd\t", $lh_bl_result{subjectend}, "\n";
			print "\n";
			print LOG "\nLeft hand best hit/hsp etc:\n";
			print LOG "\tQuery\t",$lh_bl_result{queryname}, "\n";  
			print LOG "\tHit\t", $lh_bl_result{hitname}, "\n";
			print LOG "\tScore\t", $lh_bl_result{score}, "\n";
			print LOG "\t%id\t", $lh_bl_result{id}, "\n";
			print LOG "\tMLength\t", $lh_bl_result{hsplength}, "\n";
			print LOG "\tQLength\t", $lh_bl_result{querylength}, "\n";
			print LOG "\tSLength\t", $lh_bl_result{hitlength}, "\n";
			print LOG "\tQStrand\t", $lh_bl_result{qstrand}, "\n";
			print LOG "\tHStrand\t", $lh_bl_result{hstrand}, "\n";
			print LOG "\tQStart\t", $lh_bl_result{querystart}, "\n";
			print LOG "\tQEnd\t", $lh_bl_result{queryend}, "\n";
			print LOG "\tSStart\t", $lh_bl_result{subjectstart}, "\n";
			print LOG "\tSEnd\t", $lh_bl_result{subjectend}, "\n";
			print LOG "\nAlignments:\n";
			$alIO->write_aln($lh_bl_result{aln});
			print LOG "\n";
		}
		
		
		#blast the cap contigs versus the right hand end! Store the best result with a bunch of stats...
		$cmd = "blastn -db right_hand -query forCAPall.fa.cap.contigs -evalue 0.00001 -outfmt 6";
		print  "\n\n$cmd\n";
		system($cmd);
		
		my %rh_bl_result;
		$rh_bl_result{score} = 0;
		my $rh_sio = Bio::SearchIO->new(-file=>"blastn -db right_hand -query forCAPall.fa.cap.contigs -evalue 0.00001 |", -format=>'blast');
		while (my $res = $rh_sio->next_result) {
			next if $res->no_hits_found;
			while (my $hit = $res->next_hit) {
				while (my $hsp = $hit->next_hsp) {
					if($hsp->bits > $rh_bl_result{score}){
						#we have a better hsp!  Store it!
						$rh_bl_result{score} = $hsp->bits;
						$rh_bl_result{queryname} = $res->query_name;
						$rh_bl_result{hitname} = $hit->name;
						$rh_bl_result{id} = $hsp->percent_identity;
						$rh_bl_result{hsplength} = $hsp->length;
						$rh_bl_result{querylength} = $res->query_length;
						$rh_bl_result{hitlength} = $hit->length;
						$rh_bl_result{qstrand} = $hsp->strand("query");
						$rh_bl_result{hstrand} = $hsp->strand("hit");
						$rh_bl_result{querystart} = $hsp->query->start;
						$rh_bl_result{queryend} = $hsp->query->end;
						$rh_bl_result{subjectstart} = $hsp->sbjct->start;
						$rh_bl_result{subjectend} = $hsp->sbjct->end;
						$rh_bl_result{aln} = $hsp->get_aln;
					}
				}
			}
		}
		if($rh_bl_result{score} == 0){
			print "\nNo blast hits found on left hand end!\n";
			print LOG "\nNo blast hits found on left hand end!\n";
		} else {
			print "\nRight hand best hit/hsp etc:\n";
			print "\tQuery\t",$rh_bl_result{queryname}, "\n";  
			print "\tHit\t", $rh_bl_result{hitname}, "\n";
			print "\tScore\t", $rh_bl_result{score}, "\n";
			print "\t%id\t", $rh_bl_result{id}, "\n";
			print "\tMLength\t", $rh_bl_result{hsplength}, "\n";
			print "\tQLength\t", $rh_bl_result{querylength}, "\n";
			print "\tSLength\t", $rh_bl_result{hitlength}, "\n";
			print "\tQStrand\t", $rh_bl_result{qstrand}, "\n";
			print "\tHStrand\t", $rh_bl_result{hstrand}, "\n";
			print "\tQStart\t", $rh_bl_result{querystart}, "\n";
			print "\tQEnd\t", $rh_bl_result{queryend}, "\n";
			print "\tSStart\t", $rh_bl_result{subjectstart}, "\n";
			print "\tSEnd\t", $rh_bl_result{subjectend}, "\n";
			print "\n";
			print LOG "\nRight hand best hit/hsp etc:\n";
			print LOG "\tQuery\t",$rh_bl_result{queryname}, "\n";  
			print LOG "\tHit\t", $rh_bl_result{hitname}, "\n";
			print LOG "\tScore\t", $rh_bl_result{score}, "\n";
			print LOG "\t%id\t", $rh_bl_result{id}, "\n";
			print LOG "\tMLength\t", $rh_bl_result{hsplength}, "\n";
			print LOG "\tQLength\t", $rh_bl_result{querylength}, "\n";
			print LOG "\tSLength\t", $rh_bl_result{hitlength}, "\n";
			print LOG "\tQStrand\t", $rh_bl_result{qstrand}, "\n";
			print LOG "\tHStrand\t", $rh_bl_result{hstrand}, "\n";
			print LOG "\tQStart\t", $rh_bl_result{querystart}, "\n";
			print LOG "\tQEnd\t", $rh_bl_result{queryend}, "\n";
			print LOG "\tSStart\t", $rh_bl_result{subjectstart}, "\n";
			print LOG "\tSEnd\t", $rh_bl_result{subjectend}, "\n";
			print LOG "\nAlignments:\n";
			$alIO->write_aln($rh_bl_result{aln});
			print LOG "\n";
		}
		
		#Rule tree for interpreting blast output!
		
		#If LH and RH blast results indicate same query contig then possible join!
		if($lh_bl_result{score} && $rh_bl_result{score} && $lh_bl_result{queryname} eq $rh_bl_result{queryname}){
			print "Have a possible join between $lh_bl_result{hitname} and $rh_bl_result{hitname}\n";
			print LOG "Have a possible join between $lh_bl_result{hitname} and $rh_bl_result{hitname}\n";
			#are the strands the same..
			if($lh_bl_result{hstrand} eq $rh_bl_result{hstrand}){
				if($lh_bl_result{id} >= $min_match && $lh_bl_result{hsplength} >= $min_match_length){
					if($rh_bl_result{id} >= $min_match && $rh_bl_result{hsplength} >= $min_match_length){
						print "Join meets quality requirements!\n";
						print LOG "Join meets quality requirements!\n";
						if($lh_bl_result{querystart} > $rh_bl_result{queryend} || $lh_bl_result{queryend} > $rh_bl_result{querystart}){
							print "Possibility of an overlap between left and right! No join sequence presented here, you'll need to look at this one manually..\n";
							print LOG "Possibility of an overlap between left and right! No Join sequence presented here, you'll need to look at this one manually..\n";
						}
						else {
							$joins{$lh_bl_result{hitname}}->{right} = $rh_bl_result{hitname};
							my $joinseq;
							if($lh_bl_result{hstrand} > 0){
								$joinseq = substr($queries{$lh_bl_result{queryname}}, $lh_bl_result{queryend}, ($rh_bl_result{querystart} - 1 - $lh_bl_result{queryend}));
							}
							else {
								$joinseq = &revcomp(substr($queries{$lh_bl_result{queryname}}, $lh_bl_result{queryend}, ($rh_bl_result{querystart} - 1 - $lh_bl_result{queryend})));
							}
							print "Join sequence: $joinseq\n";
							$logfile =~ m/(\d+)/;
							print LOG "Join sequence:\n>gap_${1}_join\n$joinseq\n";
							$joins{$lh_bl_result{hitname}}->{join} = $joinseq;
						}
					}
					else {
						print "Right hand side blast match failed quality checks!  Not a valid join!\n";
						print LOG "Right hand side blast match failed quality checks!  Not a valid join!\n";
					}
				}
				else {
					print "Left hand side blast match failed quality checks!  Not a valid join!\n";
					print LOG "Left hand side blast match failed quality checks!  Not a valid join!\n";
				}
			}
			else { #need to report different strands, not a joining sequence.. possible inverted repeat.
				print "Not a join.. LH on different strand to RH.. \n";
				print LOG "Not a join.. LH on different strand to RH.. \n";
			}
		}
		else { #different contigs at either end and so extensions possible (no join)!
			print "No join! Possible extensions!\n";
			print LOG "No join! Possible extensions!\n";
			#now look at left hand side first!
			#does a match exist and does it pass quality checks?
			if($lh_bl_result{score} && $lh_bl_result{id} >= $min_match && $lh_bl_result{hsplength} >= $min_match_length){
				print "left hand extension looks ok\n";
				print LOG "left hand extension looks ok\n";
				#look for extension sequence...
				if($lh_bl_result{querylength} > $lh_bl_result{queryend}){
					my $extension;
					if($lh_bl_result{hstrand} > 0){
						$extension = substr($queries{$lh_bl_result{queryname}}, $lh_bl_result{queryend});
					}
					else{
						$extension = &revcomp(substr($queries{$lh_bl_result{queryname}}, $lh_bl_result{queryend}));
					}
					print "Left hand side extension: $extension\n";
					print LOG "Left hand side extension:\n";
					$logfile =~ m/(\d+)/;
					print LOG ">gap_${1}_left\n$extension\n";
					$scaffolds{$lh_bl_result{hitname}} .= $extension;
				}
				else {
					print "No extension of left hand side possible.  No contig overhang into gap.\n";
					print LOG "No extension of left hand side possible.  No contig overhang into gap.\n";
				}
			}
			else {
				print "Left hand side extension failed quality checks.\n";
				print LOG "Left hand side extension failed quality checks.\n";
			}
			if($rh_bl_result{score} && $rh_bl_result{id} >= $min_match && $rh_bl_result{hsplength} >= $min_match_length){
				print "right hand extension looks ok\n";
				print LOG "right hand extension looks ok\n";
				#look for extension sequence...
				if($rh_bl_result{querystart} > 1){
					my $extension;
					if($rh_bl_result{hstrand} > 0){
						$extension = substr($queries{$rh_bl_result{queryname}}, 0, $rh_bl_result{querystart} - 1);
					}
					else {
						$extension = &revcomp(substr($queries{$rh_bl_result{queryname}}, 0, $rh_bl_result{querystart} - 1));
					}
					print "Right hand side extension: $extension\n";
					print LOG "Right hand side extension:\n";
					$logfile =~ m/(\d+)/;
					print LOG ">gap_${1}_right\n$extension\n";
					my $tmp = $scaffolds{$rh_bl_result{hitname}};
					$scaffolds{$rh_bl_result{hitname}} = $extension . $tmp;
				}
				else {
					print "No extension of right hand side possible.  No contig overhang into gap.\n";
					print LOG "No extension of right hand side possible.  No contig overhang into gap.\n";
				}
			}
			else {
				print "Right hand side extension failed quality checks.\n";
				print LOG "Right hand side extension failed quality checks.\n";
			}
			
		}
		print LOG "\n";
		close LOG;
		$count ++;
	}
	
	#output the new scaffolds
	my $outfile = "new_scaffolds_${iter_count}.fa";
	my $newscaffio = Bio::SeqIO->new(-file => ">$outfile", -format => 'Fasta');
	open JOINS, ">>joins.txt";
	for(my $i = 0; $i < @scaff_order; $i ++){
		if($joins{$scaff_order[$i]}){
			print "There is a join next to $scaff_order[$i]\n";
			print JOINS "There is a join next to $scaff_order[$i]\n";
			my $newid = $scaff_order[$i] ."_". $scaff_order[$i+1];
			my $newseq = $scaffolds{$scaff_order[$i]} . $joins{$scaff_order[$i]}->{join} . $scaffolds{$joins{$scaff_order[$i]}->{right}};
			my $seq = Bio::Seq->new(-id => $newid, -seq => $newseq);
			$newscaffio->write_seq($seq);
			$i ++;
		}
		else {
			my $seq = Bio::Seq->new(-id => $scaff_order[$i], -seq => $scaffolds{$scaff_order[$i]});
			$newscaffio->write_seq($seq);
		}
	}

	$iter_count ++;
	$num_iters --;
	
	print  "Num_iters at end of iteration is $num_iters\n";
	
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"i|iterations=i", VAR=>\$num_iters, DEFAULT=>1, DESC=>"The number of iterations to perform"},
		{OPT=>"d|dir=s", VAR=>\$dir, DEFAULT=>"out_dir", DESC=>"The output directory for the iterations"},
		{OPT=>"f|flank=i", VAR=>\$flank, DEFAULT=>300, DESC=>"The size of the flanking seqeunce to use"},
		{OPT=>"t|threads=i", VAR=>\$cpus, DEFAULT=>4, DESC=>"The number of threads to use for the read mapping steps"},
		{OPT=>"k|k_val=i", VAR=>\$k_val, DEFAULT=>31, DESC=>"The k-value to use in the assembly steps"},
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
	print STDERR "Usage: $0 [options] scaffolds_file.fa reads_files.fa\n";
	foreach (@Options) {
		printf STDERR "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------

#revcomp!
sub revcomp {
	my $x = shift;
	$x =~ tr/ATCGatcg/TAGCtagc/;
	my $rev = reverse $x;
	return $rev;
}

#----------------------------------------------------------------------
# Read and Write a Fasta File
#
#my $in   = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
#my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
#while (my $seq = $in->next_seq) {
#  $out->write_seq($seq);
#} 

#----------------------------------------------------------------------
# Parse a Genbank File

#my $gbk = Bio::SeqIO->new(-file=>'input.gbk', -format=>'genbank');
#while (my $seq = $gbk->next_seq) {
#  for my $f ($seq->get_SeqFeatures) {
#    next if $f->primary_tag eq 'source';
#    print $f->primary_tag;
#    my($tag) = $f->has_tag('locus_tag') ? 
#    my($tag) = $f->get_tag_values('locus_tag');      
#    $tag ||= '(no tag)';
#    print "\t", $tag; 
#    print "\t", $f->location->gff_string;
#    print "\n";
#  }
#}

#----------------------------------------------------------------------
# Parse a Blast result file

#my $bls = Bio::SearchIO->new(-file=>'out.blast', -format=>'blast');
#while (my $res = $bls->next_result) {
#  next if $res->no_hits_found;
#  print STDERR $res->query_name,"\n";
#  while (my $hit = $res->next_hit) {
#    print STDERR "\t", $hit->name,"\n";
#    while (my $hsp = $hit->next_hsp) {
#      print STDERR "\t\t", $hsp->significance,"\n";
#    }
#  }
#}

