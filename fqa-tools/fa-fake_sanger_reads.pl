#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use List::Util qw(max min);

my(@Options, $verbose, $readlen, $insert, $stride, $library);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $nread=0;
my $nwrote=0;
my $nbp;

while (my $seq = $in->next_seq) {
  $nread++;
  print STDERR "Sangerizing: #$nread ",$seq->id," ",$seq->length, " bp.\n";
  sangerize($seq, $nread);
  $nbp += $seq->length;
} 
print STDERR "Read $nread sequences ($nbp bp).\n";
print STDERR "Wrote $nwrote reads (", $readlen*$nwrote, " bp).\n";
printf STDERR "Coverage is %.2fx\n", $readlen*$nwrote/$nbp;

#----------------------------------------------------------------------

# >kbrh1_224fr_a01.ab1 template=kbrh1_224f~_a01.ab1 dir=r library=SimonA09

sub sangerize {
  my($seq, $prefix) = @_;
  my $count=0;
  my $L = $seq->length;
  for (my $start=1; $start < $L-$readlen-$insert-$readlen; $start+=$stride) {
    $count++;
    my $template = "read_${prefix}_${count}";
    print STDERR "Read-pair: $template\n" if $verbose;
    # left mate
    my $left = $seq->trunc($start, $start+$readlen-1);
    $left->id("$template.b");
    $left->desc("template=$template dir=f library=$library");
    # right mate (revcom!)
    my $right = $seq->trunc($insert+$start, $insert+$start+$readlen-1)->revcom;
    $right->id("$template.g");
    $right->desc("template=$template dir=r library=$library");
    # spit out
    $out->write_seq($left, $right);
    $nwrote+=2;
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"l|library=s",  VAR=>\$library, DEFAULT=>'Sanger', DESC=>"Library name"},
    {OPT=>"r|readlen=i",  VAR=>\$readlen, DEFAULT=>800, DESC=>"Read length"},
    {OPT=>"i|insert=i",  VAR=>\$insert, DEFAULT=>3000, DESC=>"Insert size"},
    {OPT=>"s|stride=i",  VAR=>\$stride, DEFAULT=>100, DESC=>"Stride/overlap between reads"},
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
  print "Usage: $0 [options] contigs.fasta > fake_sanger_reads.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
