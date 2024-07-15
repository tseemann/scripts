#!/usr/bin/perl -w
use strict;
use Fatal;
use Bio::SeqIO;
#use File::Spec;
use File::Temp qw(tempdir);
use Cwd qw(getcwd abs_path);
use Data::Dumper;

my $PYROBAYES_EXE = 'PyroBayes';
my $SFFINFO_EXE = 'sffinfo';
my $STEM = '454';

my(@Options, $verbose, $pair_file, $sing_file, $crap_file, $bayes, 
             $lsuffix, $rsuffix, $minlen, $auto);
setOptions();

if ($auto) {
  $pair_file = 'p.fa';
  $sing_file = 's.fa';
  $crap_file = 'c.fa';
}

($pair_file or $sing_file) or usage();

require_exe('cross_match', 1);

my $dir = getcwd();
my @sff = map { abs_path($_) } @ARGV;

my $tmpdir = tempdir(CLEANUP=>1);
chdir $tmpdir;

if ($bayes) {
  require_exe($PYROBAYES_EXE);
  for my $sff (@sff) {
    run_cmd("$PYROBAYES_EXE --outStub tmp -i '$sff'");
    run_cmd("cat tmp.fasta >> $STEM.fasta");
    run_cmd("cat tmp.fasta.qual >> $STEM.fasta.qual");
    unlink "tmp.fasta", "tmp.fasta.qual";
  }
}
else {
  require_exe('sffinfo');
  for my $sff (@sff) {
    run_cmd("sffinfo -s '$sff' >> $STEM.fasta");
    run_cmd("sffinfo -q '$sff' >> $STEM.fasta.qual");
  }
}

# http://seqanswers.com/forums/showthread.php?t=2698
my $FLX = 'GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC';
my $TIT = 'TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG';
open my $lfh, '>', "$STEM.screen";
print $lfh ">454_FLX\n$FLX\n";
print $lfh ">454_Titanium\n$TIT\n";
close $lfh;
run_cmd("cat $STEM.screen");

run_cmd("cross_match $STEM.fasta $STEM.screen -screen > screen.stdout 2> screen.stderr");

my $in = Bio::SeqIO->new(-file=>"$STEM.fasta.screen", -format=>'Fasta');

print STDERR "Pair output file: $pair_file\n";
open my $pfh, '>', "$dir/$pair_file";
my $npair=0;

print STDERR "Singleton output file: $sing_file\n";
open my $sfh, '>', "$dir/$sing_file";
my $nsing=0;

print STDERR "Crap output file: $crap_file\n";
open my $cfh, '>', "$dir/$crap_file";
my $ncrap=0;

my $nread=0;
my %parts;
my %linkers;

while (my $seq = $in->next_seq) {
  my $id = $seq->display_id;
  $nread++;

  my @s = grep { $_ } ( split m/(X+)/i, uc($seq->seq) );
  $parts{scalar(@s)}++;
  my $numlinker = scalar( grep { m/X/ } @s );
  $linkers{$numlinker}++;
  
  for my $x (@s) {
    die ">$id: split has undef in it!\n".$seq->seq if not defined $x;
    die ">$id: split has zero length sequence in it!\n".$seq->seq if $x eq '';    
  }
  
  if (@s==3 and $s[0]!~/X/ and $s[1]=~/^X+$/ and $s[2]!~/X/) {
    if (length $s[0] < $minlen and length $s[2] < $minlen) {
      # both too short, put in crapper
      print $cfh ">$id$lsuffix\n$s[0]\n";
      print $cfh ">$id$rsuffix\n$s[2]\n";
      $ncrap+=2;
    }
    elsif (length $s[0] < $minlen) {
      # left too short, right ok as single
      print $cfh ">$id$lsuffix\n$s[0]\n";
      $ncrap++;
      print $sfh ">$id$rsuffix\n$s[2]\n";
      $nsing++;
    }
    elsif (length $s[2] < $minlen) {
      # right too short, left ok as single
      print $sfh ">$id$lsuffix\n$s[0]\n";
      $nsing++;
      print $cfh ">$id$rsuffix\n$s[2]\n";
      $ncrap++;
    }
    else {
      ## we revcom both as 454 is "<= =>" orientation, velvet expects "=> <="
      #print $pfh ">$id/1\n", revcom($s[0]), "\n";
      #print $pfh ">$id/2\n", revcom($s[2]), "\n";
      # http://seqanswers.com/forums/showthread.php?t=8677

      # not sure if FLX and Titanium differ here!
      # we revcom R as 454 is "=> =>" orientation, velvet/gap4 expects "=> <="    
      #print $pfh ">$id$lsuffix\n",        $s[0] , "\n";
      #print $pfh ">$id$rsuffix\n", revcom($s[2]), "\n";

      # 5 DEC 2011
      # http://seqanswers.com/forums/showthread.php?t=8677
      print $pfh ">$id$lsuffix\n", revcom($s[0]) , "\n";
      print $pfh ">$id$rsuffix\n",        $s[2]  , "\n";

      # count them
      $npair++;
    }
  }
  elsif (@s < 3) {
    die "too many linkers?" if $numlinker > 1;
    print $sfh ">$id\n", (grep { !m/X/ } @s), "\n";
    $nsing++;
  }
  else {
    print $cfh ">$id\n", $seq->seq, "\n";
    $ncrap++;
  }
} 

#print Dumper(\%parts, \%linkers);

printf STDERR "Read: %d reads\n", $nread;
printf STDERR "Wrote: %d pairs + %d singeltons + %d crap\n", $npair, $nsing, $ncrap;

# change back to original dir so File::Temp can clean up
chdir $dir;

# END OF SCRIPT

#-------------------------------------------------------------------------

sub run_cmd {
  print STDERR "Running: @_\n";
  if (system(@_)) {
    print STDERR "ERROR $? : $!\n";
    exit $?;
  }
}

#-------------------------------------------------------------------------

sub revcom {
  my $dna = shift;
  $dna = reverse($dna);
  $dna =~ tr/ACGTacgt/TGCAtgca/;
  return $dna;
}  

#-------------------------------------------------------------------------

sub require_exe {
  my($bin, $fatal) = @_;
  for my $dir (File::Spec->path) {
    my $exe = File::Spec->catfile($dir, $bin);
    return $exe if -x $exe; 
  }
  die "Could not find '$bin' - please install it." if $fatal;
  return;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"h|help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"a|auto!",  VAR=>\$auto, DEFAULT=>0, DESC=>"Auto-enable -s s.fa -p p.fa -c c.fa"},
    {OPT=>"b|bayes!",  VAR=>\$bayes, DEFAULT=>0, DESC=>"Use PyroBayes instead of sffinfo"},
    {OPT=>"p|pairs=s",  VAR=>\$pair_file, DEFAULT=>'/dev/null', DESC=>"Output interleaved pairs"},
    {OPT=>"c|crap=s",  VAR=>\$crap_file, DEFAULT=>'/dev/null', DESC=>"Output crap/short sequences"},
    {OPT=>"s|singletons=s",  VAR=>\$sing_file, DEFAULT=>'/dev/null', DESC=>"Output singletons"},
    {OPT=>"l|lsuffix=s",  VAR=>\$lsuffix, DEFAULT=>'/1', DESC=>"Left read name suffix"},
    {OPT=>"r|rsuffix=s",  VAR=>\$rsuffix, DEFAULT=>'/2', DESC=>"Right read name suffix"},
    {OPT=>"m|minlen=i",  VAR=>\$minlen, DEFAULT=>24, DESC=>"If less than this, put in the -c crapper"},
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
  print "Usage: $0 [options] -p pairs.fasta -s singletons.fasta File1.sff [File2.sff ...]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__
>AdaptorA
CTGAGACAGGGAGGGAACAGATGGGACACGCAGGGATGAGATGG
>AdaptorB
CTGAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGG
>5prime454adaptor
GCCTCCCTCGCGCCATCAGATCGTAGGCACCTGAAA
>3prime454adaptor
GCCTTGCCAGCCCGCTCAGATTGATGGTGCCTACAG
>BanOne
CATGATTGATGGTGCCTACAG
>BanTwo
CATGATCGTAGGCACCTGAAA
>FLX_linker
GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC
>Titanium_1
TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG
>Titanium_2
CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA
>OligoLinker
CTCGAGAATTCTGGATCCTC


GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC
TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG
CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA

TTTGTTTTGTTTGTTTTGTTTT
AAAACAAAACAAACAAAACAAA
AAAAAAAAAAAAAAAAAAAA
TCCTGACGTCAAGTGATCCATT
../src/SeqTrim.cpp
CTCGAGAATTCTGGATCCTC
trimmedReads.fna
overlapsVec.at(i).size() == 1
Setting up trimming/screening overlap detection...
Reading screening database: %s
Reading trimming database: %s
        Warning: Suspected 5' primer %s, %d exact matches found.
                Warning: Suspected 3' primer %s, %d exact matches found.
                AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                AAGCAGTGGTATCAACGCAGAGTACGCGGG
                cdnaSeqCache == __null && cdnaSmallOverlapper == __null
                nimblegenVecSeqCache == __null && nimblegenVecOverlapper == __null
                TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATC
                AGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCG
                AAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAG
                CTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCC
                AGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATC
                AGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGA
                CTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCG
                TTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGG
                CGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTT
                GTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATT
                AAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGA
                GGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGG
                GAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAA
                AAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTG
                AGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATG
                TAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAA
                GCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTAT
                CACGAGGCCCTTTCGTC
                pairedVecSeqCache == __null && pairedVecOverlapper == __null
                CTGAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGG
                



