#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::Index::Fasta;

my(@Options, $verbose, $flank, $dir, $contigfn, $maskfn, $pairfn, $lower, $upper);
setOptions();
$|=1;

$dir or die "need --dir <outdir>";
-d $dir and die "whoops, --dir $dir already exists!";
mkdir $dir or die "crap, i can't mkdir $dir :-(";

-r $contigfn or die "can't read contigs '$contigfn'";
-r $maskfn or die "can't read screen/mask sequences '$maskfn'";
-r $pairfn   or die "can't read mate-pair reads '$pairfn'";
$lower or die "please supply --lower";
$upper or die "please supply --upper";
$upper > $flank and die "max insert size bigger than flank - you sure?";
$upper > $lower or die "upper needs to be bigger than lower";

my $blastdb = "$dir/blast.db";
print STDERR "Making blastn database: $blastdb\n";
qx(formatdb -i $contigfn -p F -n $blastdb -l /dev/null);

my $idxfn = "$dir/pairs.idx";
print STDERR "Indexing '$pairfn' to '$idxfn'. Please wait....\n";
my $idx = Bio::Index::Fasta->new(-filename=>$idxfn, -write_flag=>1);
$idx->make_index($pairfn);

my $pairdb = "$dir/pairs.db";
print STDERR "Formatdb '$pairfn' to '$pairdb'. Please wait....\n";   
qx(formatdb -i $pairfn -p F -n $pairdb -l /dev/null);

print STDERR "Screening contigs in '$contigfn' from $maskfn\n";
qx(cross_match $contigfn $maskfn -screen 2> /dev/null);
qx(mv $contigfn.screen $dir/);
my $screenfn = "$dir/$contigfn.screen";

my $webfn = "$dir/index.html";
print STDERR "Creating: $webfn\n";
open my $webfh, '>', $webfn;
my($pwd) = qx(pwd); chomp $pwd;
print $webfh "<h1>$contigfn</h1>\n";

my $cin = Bio::SeqIO->new(-file=>$contigfn, -format=>'Fasta');
my $min = Bio::SeqIO->new(-file=>$screenfn, -format=>'Fasta');
while (my $cseq = $cin->next_seq) {
  my $mseq = $min->next_seq;
  my $id = $cseq->display_id;
#  my $cs = uc $cseq->revcom->seq; # WHY AM I REVCOMMING ????
  my $cs = uc $cseq->seq; # WHY AM I REVCOMMING ????
  if (1) {
    my $s = $mseq->seq;
    $s =~ s/X/N/gi;
    $mseq->seq($s); 
  }
#  my $ms = uc $mseq->revcom->seq; # WHY REVCOM???
  my $ms = uc $mseq->seq; # WHY REVCOM???
  print $webfh "<hr>\n";
  find_mates( "$id.Left" ,substr($cs,0,$flank), substr($ms,0,$flank), $webfh );
  find_mates( "$id.Right", substr($cs,-$flank,$flank), substr($ms,-$flank,$flank), $webfh );
} 

chdir $dir && qx(ln -s . $dir);

#----------------------------------------------------------------------

sub find_mates {
  my($id, $cseq, $mseq, $webfh) = @_;

  my $name = $id;
  $name =~ s/[^\w0-9_.-]/_/g;
  
  my $prefix = "[$name]";

  my $sfn = "$dir/$name.fasta";
  print STDERR "$prefix Writing flank: $sfn\n";
  my $sout = Bio::SeqIO->new(-file=>">$sfn", -format=>'fasta');
  $sout->write_seq( Bio::Seq->new(-id=>$id, -seq=>$cseq) );
  undef $sout;

  my $msfn = "$dir/$name.fasta.screen";
  print STDERR "$prefix Writing screened flank: $msfn\n";
  my $mout = Bio::SeqIO->new(-file=>">$msfn", -format=>'fasta');
  $mout->write_seq( Bio::Seq->new(-id=>"$id.screen", -seq=>$mseq) );
  undef $mout;
  
  print STDERR "$prefix BLASTing for mates ... ";
  my @ids = qx(blastall -a 4 -v 9999 -b 9999 -p blastn -i $msfn -e 1E-11 -m 8 -d $pairdb | cut -f2 | cut -d/ -f1 | sort | uniq);
  chomp @ids;
  my $npairs = scalar @ids;
  print STDERR "found $npairs pairs.\n";

  my $rfn = "$dir/$name.reads";
  qx(cat $sfn > $rfn); # put original in read set

  my $cfn = "$rfn.con";  # constraints!!!
  open my $cfh, '>', $cfn;

  my $reads=1;
  my $rfh = Bio::SeqIO->new(-file=>">>$rfn", -format=>'fasta');
  print STDERR "$prefix Collecting mates: $rfn\n";
  for my $id (@ids) {
    for my $p (1 .. 2) {
      my $s = $idx->fetch("$id/$p") or next;
      my $cap3id = "$id.".($p==1 ? 'f' : 'r');
      $s->display_id($cap3id);
      $rfh->write_seq($s);
      print $cfh "$cap3id\t";
      $reads++;
    }
    print $cfh "$lower\t$upper\n";
  }
  close $cfh;
  
  print STDERR "$prefix Assembling $reads reads in: $rfn\n";
  qx(cap3 $rfn);
  my($contigs) = qx(grep -c '>' $rfn.cap.contigs); chomp $contigs;
  my($singlets) = qx(grep -c '>' $rfn.cap.singlets); chomp $singlets;
  print STDERR "Assembled into $contigs contigs and $singlets singlets.\n";
  
  print STDERR "$prefix BLAST contigs against assembly\n";
  qx(blastall -T T -p blastn -i $rfn.cap.contigs -o $rfn.cap.contigs.bls -d $blastdb -e 1E-11 -a 4);
  
  print STDERR "$prefix Cleaning up\n";
  unlink "$rfn.cap.contigs.links", "$rfn.cap.contigs.qual", 
         "$rfn.cap.info", "$rfn.cap.ace";
  
  $pwd='.';
  print $webfh <<"EOFWEB";
<h2>$id</h2>
<ul>
<li><a href='$sfn'>flanking sequence</a>
<li><a href='$msfn'>masked flanking sequence</a>
<li><a href='$rfn'>anchored mate pairs</a> ($npairs)
<li><a href='$cfn'>mate pair constraints</a> ($npairs)
<li><a href='$rfn.cap.contigs'>contigs</a> ($contigs)
<li><a href='$rfn.cap.singlets'>singlets</a> ($singlets)
<li><a href='$rfn.cap.contigs.bls'>blastn: contigs vs assembly</a>
</ul>  

EOFWEB
  
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"d|dir=s",  VAR=>\$dir, DEFAULT=>'', DESC=>"Output directory"},
    {OPT=>"c|contigs=s",  VAR=>\$contigfn, DEFAULT=>'', DESC=>"Ordered contigs to find links between [FASTA]"},
    {OPT=>"m|maskfile=s",  VAR=>\$maskfn, DEFAULT=>'', DESC=>"Repeats sequences to mask [FASTA]"},
    {OPT=>"p|pairs=s",  VAR=>\$pairfn, DEFAULT=>'', DESC=>"Interleaved mate-pair reads [FASTA]"},
    {OPT=>"l|lower=i",  VAR=>\$lower, DEFAULT=>0, DESC=>"Minimum insert size"},
    {OPT=>"u|upper=i",  VAR=>\$upper, DEFAULT=>0, DESC=>"Maximum insert size"},
    {OPT=>"f|flank=i",  VAR=>\$flank, DEFAULT=>3000, DESC=>"Flank size"},
  );

  (@ARGV) or (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
