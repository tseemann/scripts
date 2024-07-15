#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use Bio::SearchIO;

my(@Options, $verbose, $folder, $contig_file, $read_file, 
             $dist, $len, $minsize, $extend);
setOptions();

#$folder or die "Please specify a destination folder with --folder";
#mkdir $folder or die "Unable to create folder '$folder': $!";
-s $contig_file or die "Contigs file '$contig_file' is empty, use --contigs <fasta>";
-s $read_file or die "Reads file '$read_file' is empty, use --reads <fasta>";
$minsize ||= ($dist + $len); # need at least this many bases!

my $in = Bio::SeqIO->new(-file=>$contig_file, -format=>'Fasta');
my $nctg=0;
while (my $seq = $in->next_seq) {
  next unless $seq->length >= $minsize;
  process_contig($seq);
  $seq = $seq->revcom;
  $seq->id($seq->id . "c");
  process_contig($seq);
  $nctg++;
} 
print "Processed $nctg contigs >= $minsize bp.\n";

#......................................................................

sub process_contig {
  my($seq) = @_;
  my $id = $seq->id;
  my $L = $seq->length;
  printf STDERR "Processing: %s | %d bp.\n", $id, $L if $verbose;
  my $oligo = uc substr $seq->seq, -($dist+$len), $len;
  print '-'x70, "\n";
  my $left = $L - $dist - $len;
  print "$id | $L bp | 5'-($left)-$oligo-($dist)-3'\n";
  my @read = qx(grep -F -i --mmap '$oligo' $read_file);
  chomp @read;
  @read = map { uc } @read;
  map { s/(^.*?$oligo)// } @read;
  if ($extend) {
    # only show reads that extend the contig when -x enabled
    @read = grep { length($_) > $dist } @read;
  }
  map { s/(^.{1,$dist})/\L$1/ } @read;
  @read = grep { $_ } @read; # remove empty read strings
  @read = reverse sort @read;
  print map { "$_\n" } @read;
  print "(no reads matched)\n" unless @read;
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"this help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"verbose output"},
#    {OPT=>"f|folder=s",  VAR=>\$folder, DEFAULT=>'', DESC=>"folder to put everything"},
    {OPT=>"c|contigs=s",  VAR=>\$contig_file, DEFAULT=>'/dev/null', DESC=>"input contig (.fasta)"},
    {OPT=>"r|reads=s",  VAR=>\$read_file, DEFAULT=>'/dev/null', DESC=>"input reads file (.fasta/.fastq/.raw)"},
    {OPT=>"d|dist=i",  VAR=>\$dist, DEFAULT=>10, DESC=>"distance from end of contig"},
    {OPT=>"l|len=i",  VAR=>\$len, DEFAULT=>15, DESC=>"length of segment to match"},
    {OPT=>"m|minsize=i",  VAR=>\$minsize, DEFAULT=>0, DESC=>"minimum contig size to process"},
    {OPT=>"x|extend!",  VAR=>\$extend, DEFAULT=>0, DESC=>"only show reads that would extend the contig"},
  );

  #(!@ARGV) && (usage());

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
