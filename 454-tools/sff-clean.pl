#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use File::Spec;
use File::Temp qw(tempdir);

my(@Options, $debug);
setOptions();

require_exe('sffinfo', 1);
my $dir = tempdir(CLEANUP=>1);
my $tmpfname = "$dir/bad.ids";
print STDERR "Temporary: $tmpfname\n";

for my $sff (@ARGV) { 
  $sff = File::Spec->rel2abs($sff);
  print STDERR "Processing: $sff\n";
  open my $badfh, '>', $tmpfname;
  my $in = Bio::SeqIO->new(-file=>"sffinfo -seq '$sff' |", -format=>'fasta');
  my $bad=0;
  my $num=0;
  while (my $seq = $in->next_seq) {
    $num++;
    if ($seq->seq !~ m/^[AGTC]+$/i) {
      print $badfh $seq->id."\n";
      $bad++;
    }  
  }
  print STDERR "Read: $num\n", "Good: ",$num-$bad, "\n", "Bad: $bad\n";
  my($vol,$path,$file) = File::Spec->splitpath($sff);
  my $target = File::Spec->catpath($vol, $path, "clean.$file");
  print STDERR "Writing clean version: $target\n";
  system("sfffile -e $tmpfname -o $target $sff");
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
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
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
  print "Usage: $0 [options] file1.sff [file2.sff ...]\n";
  print "Output: clean.file1.sff [ ... ]\n";
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
