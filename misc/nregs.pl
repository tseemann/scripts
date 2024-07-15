#!/usr/bin/env perl

use strict;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;
use File::Temp qw(tempfile tempdir);
use Data::Dumper;

my (@Options,@nregs,$results,$scaffs,$contigs,$verbose,$flank,$id,$mlen,$hdir,$blastextras,$blastndatabase,$blastxdatabase);
setOptions();

$scaffs = $ARGV[0];
$contigs = $ARGV[1];

my %contigs_seqs;

unless($scaffs && -r $scaffs){ print STDERR "\nNo scaffs file specified or not readable.\n\n"; &usage; }
unless($contigs && -r $contigs){ print STDERR "\nNo contigs file specified or not readable.\n\n"; &usage; }

print STDERR "Finding N Regions in $scaffs\n" if $verbose;

@nregs = getNregions($scaffs, $flank, "N");

#foreach (@nregs){
#	print $_ . "\n";
#}

#load the contigs into a hash!

my $cio = Bio::SeqIO->new(-file=>$contigs, -format=>'Fasta');
while(my $seq = $cio->next_seq){
	$contigs_seqs{$seq->id} = $seq->seq;
}
#foreach my $key(sort keys %contigs_seqs){
#	print $key . "\n" . $contigs_seqs{$key} . "\n";
#}

#now put all the left and right regions into a file...  Torst will no doubt make me use a tmp file later on...
#to do: what if either flank sequence contains Ns?

print STDERR "Making flanking sequence fasta file\n" if $verbose;

my $fasta_fh = File::Temp->new(SUFFIX=>'.fasta.');

my $tmp_out = Bio::SeqIO->new(-file => ">" . $fasta_fh, -format=> 'Fasta');

foreach my $line (@nregs){
	my @tmp = split m/\t/, $line;
	my $left_id = $tmp[0] . "_$tmp[1]_$tmp[3]_L";
	my $right_id =  $tmp[0] . "_$tmp[1]_$tmp[3]_R";
	my $Ls = Bio::Seq->new(-id=>$left_id, -seq=>$tmp[2]);
	my $Rs = Bio::Seq->new(-id=>$right_id, -seq=>$tmp[5]);
	$tmp_out->write_seq($Ls);
	$tmp_out->write_seq($Rs);
}

#my $fn = $fasta_fh->filename;
print STDERR "Filename: ". $fasta_fh . "\n" if $verbose;

print STDERR "Making blast database from $contigs\n" if $verbose;

#make the blastdb from the contigs!
my $blastdir = tempdir(CLEANUP=>1);

system("formatdb -o T -p F -i $contigs -n $blastdir/contigs");

print STDERR "Running blastn of flanking sequences versus contigs\n" if $verbose;

my $search_results = Bio::SearchIO->new(-file => "blastn -db $blastdir/contigs -query $fasta_fh -evalue 0.00001|", -format => 'blast');

#print Dumper($search_results);

my %hits;


print "Scaffold\tContig\tContigLength\tScaffStrand\tContigStrand\tMatchLength\tContigMatchLength\tScaffHitLength\tStartContig\tEndContig\tPercentIdentity\n" if $verbose;
while(my $result = $search_results->next_result){
	## $result is a Bio::Search::Result::ResultI compliant object
	while( my $hit = $result->next_hit ) {
		## $hit is a Bio::Search::Hit::HitI compliant object
		while( my $hsp = $hit->next_hsp ) {
		## $hsp is a Bio::Search::HSP::HSPI compliant object
			if( $hsp->length('total') > $mlen ) {
				if ( $hsp->percent_identity >= $id ) {
					my $out = $hit->name.
						"\t". $hit->length.
						"\t". $hsp->strand('hit').
						"\t". $hsp->strand('query').
						"\t". $hsp->length('total').
						"\t". $hsp->length('hit').
						"\t". $hsp->length('query').
						"\t". $hsp->start('hit').
						"\t". $hsp->end('hit').
						"\t". $hsp->percent_identity;
					
					#print out hit if verbose...
					print $result->query_name, "\t$out\n" if $verbose;
										
					#load hit into data structure...
					push (@{$hits{$result->query_name}}, $out);
						
				}
			}
		}  
	}
}

print STDERR "Searching hits for possible matches and producing output..\n" if $verbose;

my %seen;

my $contig_sections_file = File::Temp->new();

foreach my $key (sort keys %hits){
	next if $seen{$key};
	#count the number of hits in the array and if multiple then treat differently...
	my $numhits = scalar @{$hits{$key}};
	if($numhits > 1){
		print "$key: Multiple hits\n";
		foreach (@{$hits{$key}}){
			print $_ . "\n";
		}
		print "\n\n";
	}
	else {
		#print "$key\t${$hits{$key}}[0]\n";
		#if left, is there a right?
		if($key =~ /L$/){
			$seen{$key} ++;
			my $tmp = $key;
			$tmp =~ s/L$/R/;
			if($hits{$tmp}){
				#there is a right!
				$seen{$tmp} ++;
				#does the left and right have the same contig and the same strand?
				my @l = split /\t/, ${$hits{$key}}[0];
				my @r = split /\t/, ${$hits{$tmp}}[0];
				if($l[0] eq $r[0] && $l[2] == $r[2]){
					#they are the same
					print "Possible contig covering N region: $key, $tmp, $l[0]\n";
					my @x = split /_/, $key;
					print "Scaffold: $x[0]\n";
					print "N-region start: $x[1]\tLength of N region: $x[2]\n";
					print "\tContig\tContigLength\tScaffStrand\tContigStrand\tMatchLength\tContigMatchLength\tScaffHitLength\tStartContig\tEndContig\tPercentIdentity\n";
					print "Left hit: ", ${$hits{$key}}[0], "\n";
					print "Right hit: ", ${$hits{$tmp}}[0], "\n";
					my $ss = $x[1] - $l[6];
					my $se = $x[1] + $x[2] + $r[6];
					my $sl = $se - $ss;
					print "Scaffd position: Start: $ss\tEnd: $se\tLength: $sl\n";
					print "Contig position: Start: ";
					my $start;
					$start = $l[7] if ($l[2] == 1);
					$start = $r[7] if ($l[2] == -1);
					print $start;
					print "\tEnd: ";
					my $end;
					$end = $r[8] if ($l[2] == 1);
					$end = $l[8] if ($l[2] == -1);
					print $end, "\tLength: ", $end - $start;
					print "\n";
					#now if the start coordinate is higher than the end one, we have an overlap and is a special case!
					if($start > $end){
						#overlap!
						print "According to Blast output, these two Scaffold flanking sequences overlap by ";
						print (($start - $end) . " bases\n\n\n");
					}
					else {
						#not overlapping
						print "Contig sequence of area:\n";
						print "Contig is reverse complemented\n" if ($l[2] == -1);
						my $name = $l[0];
						$name =~ s/^lcl\|//;
						my $contigsection = "";
						$contigsection = substr($contigs_seqs{$name}, $start -1, $end-$start) if($l[2] == 1);
						$contigsection = &revcomp(substr($contigs_seqs{$name}, $start -1, $end-$start)) if($l[2] == -1);
						print $contigsection . "\n\n\n";
						print $contig_sections_file ">${key}_$l[0]\n";
						print $contig_sections_file "$contigsection \n";
					}
				}
			}
		}
	}
}

if($blastextras){
	print STDERR "Running blastn on all contig sections that matched versus nt.\n";
	my $blastn = &blastncontig($contig_sections_file);
	open HTML1, ">${hdir}_n" or die "Bugger - couldn't open ${hdir}_n for writing.\n$!";
	print HTML1 $blastn;
	
	print STDERR "Running blastx on all contig sections that matched versus nr.\n";
	my $blastr = &blastxcontig($contig_sections_file);
	open HTML2, ">${hdir}_r" or die "Bugger - couldn't open ${hdir}_r for writing.\n$!";
	print HTML2 $blastr;
}

print STDERR "Finished!\n" if $verbose;

#----------------------------------------------------------------------
#Subroutines

#getNregions - Original copyright Torsten Seemann - Monash University, adapted by Simon Gladman 2010.
sub getNregions {
	my $file = shift;
	my $fl = shift;
	my $chars = shift;
	my @out;
	
	my $in = Bio::SeqIO->new(-file=>$file, -format=>'Fasta');
		
	while (my $seq = $in->next_seq) 
	{
	  my $s = $seq->seq;
	
	  # pos($s) is the position (base-0) of the first non-char AFTER run (match)
	  # but remember DNA is (base-1) ... sigh
	  
	  while ($s =~ m/([N]+)/gi) {
	      push @out, join("\t", 
	      $seq->id, 
	      pos($s)-length($1)+1,  
	      substr($s, pos($s)-length($1)-$fl, $fl),
	      length($1),
	      $chars,
	      substr($s, pos($s), $fl),
	    );
	  }
	
	}
	return @out;	
}

#sub reverse complement a sequence string...
sub revcomp {
	my $x = shift;
	$x =~ tr/[ACGTacgt]/[TGCAtgca]/;
	my $out = reverse $x;
	return $out;
}

sub blastncontig {
	my $file = shift;
	my $search_results = `blastn -db $blastndatabase -num_descriptions 10 -num_alignments 2 -num_threads 8 -query $file -evalue 0.00001 -html`;
	return $search_results;
}

sub blastxcontig {
	my $file = shift;
	my $search_results = `blastx -db $blastxdatabase -num_descriptions 10 -num_alignments 2 -num_threads 8 -query $file -evalue 0.00001 -html`;
	return $search_results;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  #matchlen & glank recommL 250 and 300
  
  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
#    {OPT=>"m|matchlength=i",  VAR=>\$mlen, DEFAULT=>50, DESC=>"Minimum match length to score a match."},
#    {OPT=>"f|flank=i",  VAR=>\$flank, DEFAULT=>100, DESC=>"Number of bases on either side of N run to use in blat search."},
    {OPT=>"m|matchlength=i",  VAR=>\$mlen, DEFAULT=>250, DESC=>"Minimum match length to score a match."},
    {OPT=>"f|flank=i",  VAR=>\$flank, DEFAULT=>300, DESC=>"Number of bases on either side of N run to use in blat search."},
    {OPT=>"i|identity=f",  VAR=>\$id, DEFAULT=>95, DESC=>"Minimum identity to score a match."},
    {OPT=>"h|htmloutfile=s",  VAR=>\$hdir, DEFAULT=>'blastresults.html', DESC=>"HTML output with file name. Useless without --blastextras (-b)."},
    {OPT=>"b|blastextras!",  VAR=>\$blastextras, DEFAULT=>0, DESC=>"Run extra blasts of matches versus nr and nt."},
    {OPT=>"n|blastn_database=s",  VAR=>\$blastndatabase, DEFAULT=>"nt", DESC=>"database to blastn contig sections against."},
    {OPT=>"r|blastr_database=s",  VAR=>\$blastxdatabase, DEFAULT=>"nr", DESC=>"database to blastx contig sections against."},
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print STDERR "Usage: $0 [options] <Scaffolds_file.fasta> <Contigs_file.fasta> > New_scaffolds_file.fasta\n";
  foreach (@Options) {
    printf STDERR "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
