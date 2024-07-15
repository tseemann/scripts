#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Temp qw(tempdir);

my(@Options, $verbose, $ref, $input, $append, $blastn, $evalue);
setOptions();

-r $ref or die "Please supply a reference FASTA file with --ref / -r";
-r $input or die "Please supply an input FASTA file with --input / -i";

my $tmpdir = tempdir(CLEANUP => 1);
my $cmd = "formatdb -i '$ref' -p F -n $tmpdir/reference -l /dev/null";
print STDERR "Running: $cmd\n";
system($cmd);

my %flip;

$cmd = $blastn ? "blastall -p blastn" : "megablast";
$cmd .= " -i '$input' -d $tmpdir/reference -m 8 -v 1 -b 1 -e $evalue";
print STDERR "Running: $cmd\n"; 
my $bls = Bio::SearchIO->new(-file=>"$cmd |", -format=>'blasttable') or die $!;

while (my $res = $bls->next_result) {
HIT:
  while (my $hit = $res->next_hit) {
     while (my $hsp = $hit->next_hsp) {
        die "query strand should be +1 always!" if $hsp->strand('query') != 1;
        # next unless $hsp->significance < 1E-6; # done in $cmd above with -e option
        my $strand = $hsp->strand('hit');
        $flip{ $res->query_name } = ( $hsp->strand('hit') != $hsp->strand('query') );
        last HIT; # take 1st decent hit we get
    }
  }  
}    
                                                                                            
my $in = Bio::SeqIO->new(-file=>$input, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $flipped = 0;
my $nread = 0;

while (my $seq = $in->next_seq) {
  if ($flip{$seq->display_id}) {
    $flipped++;
    $seq = $seq->revcom;
    $seq->display_id( $seq->display_id . $append ) if $append;
  }
  $out->write_seq($seq);
  $nread++;
} 

print STDERR "Re-oriented $flipped of $nread sequences.";
print STDERR " Appended '$append' to IDs." if $append;
print STDERR "\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",       VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!", VAR=>\$verbose, DEFAULT=>0,  DESC=>"Verbose"},
    {OPT=>"r|ref=s",    VAR=>\$ref,     DEFAULT=>'', DESC=>"Reference sequence(s) [FASTA]"},
    {OPT=>"i|input=s",  VAR=>\$input,   DEFAULT=>'', DESC=>"Input file [FASTA]"},
    {OPT=>"a|append=s", VAR=>\$append,  DEFAULT=>'', DESC=>"Append this to ID if re-oriented"},
    {OPT=>"b|blastn!", VAR=>\$blastn,  DEFAULT=>0, DESC=>"Use BLASTN instead of MEGA-BLAST"},
    {OPT=>"e|evalue=f", VAR=>\$evalue,  DEFAULT=>1E-6, DESC=>"BLAST cutoff"},
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
  print "Usage: $0 [options] --ref reference.fa --input messy.fa [--append c] > reoriented.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
