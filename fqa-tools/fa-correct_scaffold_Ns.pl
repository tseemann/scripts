#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use IO::String;
use Data::Dumper;
use Bio::SearchIO;
use File::Temp;

my(@Options, $verbose, $maxgap, $flank, $identity);
setOptions();

my $ref = shift @ARGV;
print STDERR "Reference: $ref\n";
die "Reference '$ref' unreadable or empty" unless -r $ref and -s $ref;

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $fasta = File::Temp->new(SUFFIX=>'.fasta');
print STDERR "Writing: $fasta\n";

my %inseq;

while (my $seq = $in->next_seq) {
  my $id = $seq->display_id;
  my $s = uc $seq->seq;
  print STDERR "Loaded: $id\n" if $verbose;
  $inseq{$id} = [ split m/(N+)/, $s ];  
  for (my $ni=1; $ni < @{$inseq{$id}} ; $ni+=2) {
    die unless $inseq{$id}[$ni] =~ m/^N+$/;
    print $fasta ">$id~$ni~L\n", substr($inseq{$id}[$ni-1], -$flank, $flank), "\n"; 
    print $fasta ">$id~$ni~R\n", substr($inseq{$id}[$ni+1],  0, $flank), "\n"; 
  }
} 
#system("cat $fasta");
system("cp -f $fasta ends.fasta");
#print STDERR Dumper(\%inseq);

my $extra = '';
#$extra = "-t=dnax -q=dnax";
#my $cmd = "blat \Q$ref\E \Q$fasta\E $extra -fine -noTrimA -out=blast -minIdentity=$identity /dev/stdout";
#my $cmd = "exonerate --model affine:bestfit \Q$fasta\E \Q$ref\E";
my $cmd = "exonerate --percent $identity --fsmmemory 1024".
          " --model ungapped --bestn 1 --showcigar 1 ".
	  " --query \Q$fasta\E --target \Q$ref\E";
print STDERR "Running: $cmd\n";

#system($cmd); exit;
my $bls = Bio::SearchIO->new(-file=>"$cmd |", -format=>'exonerate', -vulgar=>1);

my %pair;

select STDERR;

while (my $res = $bls->next_result) {
#  print "Query: ", $res->query_name,"\n";
  next if $res->hits <= 0;
  my($srcid, $index, $side) = split m/~/, $res->query_name;
  my $id = "$srcid.$index";
  my $hit = $res->next_hit;
  my $hsp = $hit->next_hsp;

  print $hit->name," ", $res->query_length, " ", $hit->length_aln,
        " $side ", $hsp->strand('qry'), ' ',
	$hsp->start('hit'), "..",$hsp->end('hit'), " $id\n";
  
#  next unless $res->query_length == $hit->length_aln;
  next unless $flank == $hit->length_aln;
  
  $pair{$id}{$side} = {
    NAME    => $hit->name,
    STRAND  => $hsp->strand('hit'), # qry for blat
    START   => $hsp->start('hit'),
    END     => $hsp->end('hit'),
    INDEX   => $index,
  };
  
  my $L = $pair{$id}{L};
  my $R = $pair{$id}{R};
  
  if ($L and $R and $R->{NAME} eq $L->{NAME} and $R->{STRAND}==$L->{STRAND}) {
    print "MATCH {\n";
    for my $h ($L, $R) {
      print "\t",$h->{NAME}, " ", $h->{STRAND}, " | ",
        " ",$h->{START}, "..", $h->{END}, " $id\n";  
    }
#    print "qry strand: ", $h->strand('qry'),"\n";
    my $gap = $L->{STRAND} >= 0
	    ? $R->{START} - $L->{END} - 1 
            : $L->{START} - $R->{END} - 1
	    ;
    my $old_gap = length( $inseq{$srcid}[$index] );
    printf "\tgap was %d Ns, mapping suggests %d\n", $old_gap, $gap;
    if (abs($gap) > $maxgap) {
      print "Skipping this gap as $gap too big ($maxgap)\n";
      next;
    }
    if ($gap <= 0) {
      $gap+=3 while $gap <= 0;
      print "\tmodifying non-positive gap to positive... now $gap\n";
    }
    print "}\n";
    $inseq{$srcid}[$index] = 'N'x$gap;
  }
}

#exit;
my $fout = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'fasta');
for my $id (sort keys %inseq) {
  $fout->write_seq( 
    Bio::Seq->new(-id=>$id, -seq=>join('', @{$inseq{$id}}) )
  );
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"m|maxgap=i",  VAR=>\$maxgap, DEFAULT=>7000, DESC=>"Maximum gap to change to"},
    {OPT=>"f|flank=i",  VAR=>\$flank, DEFAULT=>48, DESC=>"Flanking length to consider"},
    {OPT=>"i|identity=f",  VAR=>\$identity, DEFAULT=>65, DESC=>"Percent identity"},
  );

  (@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] reference.fna scaffolds.fna [contigs.fna ...]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
