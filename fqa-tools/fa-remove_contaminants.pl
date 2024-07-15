#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Temp qw(tempdir);

my(@Options, $verbose, $input, $contam, $output, $dirty, $blastn, $evalue, $mincov, $minlen);
setOptions();

($mincov >= 0 or $mincov <= 1) or die "--mincov must be between 0 and 1 (a proportion)";
-r $contam or die "Please supply a reference FASTA file with --contaminant";
-r $input or die "Please supply an input FASTA file with --input";

my $tmpdir = tempdir(CLEANUP => 1);
my $contamdb = "$tmpdir/contam";
#my $cmd = "formatdb -i \Q$contam\E -p F -n $contamdb -l /dev/null";
my $cmd = "makeblastdb -in \Q$contam\E -dbtype nucl -out $contamdb -logfile /dev/null";
print STDERR "Running: $cmd\n";
system($cmd);

#$cmd = $blastn ? "blastall -a 8 -p blastn" : "megablast";
#$cmd .= " -i \Q$input\E -d $contamdb -v 500 -b 500 -e $evalue -F F";
$cmd = "blastn -query \Q$input\E -db $contamdb -evalue $evalue -num_threads 8"; 
print STDERR "Running: $cmd\n"; 
print STDERR "Parsing result...\n";
my $bls = Bio::SearchIO->new(-file=>"$cmd |", -format=>'blast') or die $!;

my %is_dirty;
my $ndirty = 0;
while (my $res = $bls->next_result) {
  while (my $hit = $res->next_hit) {
     if ($hit->frac_aligned_query > $mincov) {
       $is_dirty{ $res->query_name }++ ;
       print STDERR  ++$ndirty, " ", $res->query_name, ' ', # "#", $is_dirty{$res->query_name}, " ",
         substr($hit->description,0,15)," ", $hit->length_aln, "/", $res->query_length, "\n";
       last;
     }
     else {
#       print STDERR $res->query_name, " OK\n";
     }
   }  
}    
print STDERR "\nWriting results...\n";
                                                                                            
my $in = Bio::SeqIO->new(-file=>$input, -format=>'Fasta');
my $out = Bio::SeqIO->new(-file=>">$output", -format=>'Fasta');
my $crap = Bio::SeqIO->new(-file=>">$dirty", -format=>'Fasta');

my $read=0;
my $readbp=0;
my $wrote=0;
my $wrotebp=0;
while (my $seq = $in->next_seq) {
  $read++;
  $readbp+=$seq->length;
  if ($is_dirty{$seq->id} or $seq->length < $minlen) {
    $crap->write_seq($seq);
  }
  else {
    $out->write_seq($seq);
    $wrote++;
    $wrotebp+=$seq->length;
  }
} 

printf STDERR "Read %d sequences, %d bp, %6.2f%%\n",
  $read, $readbp, (100*$readbp/$readbp) ;
printf STDERR "Wrote %d sequences, %d bp, %6.2f%%\n",
  $wrote, $wrotebp, (100*$wrotebp/$readbp) ;
printf STDERR "Removed %d sequences, %d bp, %6.2f%%\n",
  $read-$wrote, $readbp-$wrotebp, (100*($readbp-$wrotebp)/$readbp) ;

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",       VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!", VAR=>\$verbose, DEFAULT=>0,  DESC=>"Verbose"},
    {OPT=>"contam=s",    VAR=>\$contam,     DEFAULT=>'', DESC=>"Contaminant sequences [FASTA]"},
    {OPT=>"input=s",  VAR=>\$input,   DEFAULT=>'', DESC=>"Input file [FASTA]"},
    {OPT=>"output=s",  VAR=>\$output,   DEFAULT=>'/dev/stdout', DESC=>"Kept sequences [FASTA]"},
    {OPT=>"dirty=s",  VAR=>\$dirty,   DEFAULT=>'/dev/null', DESC=>"Removed sequences [FASTA]"},
    {OPT=>"blastn!", VAR=>\$blastn,  DEFAULT=>0, DESC=>"Use BLASTN instead of MEGA-BLAST"},
    {OPT=>"evalue=f", VAR=>\$evalue,  DEFAULT=>1E-6, DESC=>"BLAST cutoff"},
    {OPT=>"mincov=f",  VAR=>\$mincov,   DEFAULT=>0.1, DESC=>"Minimum proportion contam to consider bad"},
    {OPT=>"minlen=i",  VAR=>\$minlen,   DEFAULT=>0, DESC=>"Minimum contig to bother with"},
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
  print "Usage: $0 [options] -c crap.fna -i contigs.fna > clean.fna\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
