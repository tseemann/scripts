#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use List::Util qw(sum);

my(@Options, $verbose, $noid, $include);
setOptions();

my $gbf = shift @ARGV;
print STDERR "Opening Genbank reference: $gbf\n";
my $ref = Bio::SeqIO->new(-file=>$gbf, -format=>'genbank');

my $rpf = shift @ARGV;
my %snp;
my $num=0;
print STDERR "Loading SNP report: $rpf\n";
open REPORT, $rpf;
while (<REPORT>) {
  my @x = split m/\t/;
  next unless $x[1] =~ m/^\d+$/;
  #Sequence        Position in reference   Change type     Old     New     Evidence
  #C1_prokka       45655   substitution    A       G       Gx88
  my $seqid = $noid ? 'UNKNOWN' : $x[0];
  $snp{$seqid}{$x[1]}++;
  $num++;
}
print STDERR "Found ",scalar(keys %snp)," sequences in report.\n";
print STDERR "Found $num SNPs in report.\n";
  
# $include = $include ? "$include,source" : 'source';
print STDERR "Counting SNPS in: $include\n";
my $incre = join '|', (split m/\s*,\s*/, $include) ;
print STDERR "Regexp = $incre\n" if $verbose;

print STDERR "Reading Genbank sequences...\n";
my %count;
my $used=0;
my $totbp=0;
while (my $seq = $ref->next_seq) 
{
  $totbp += $seq->length;
  my $seqid = $noid ? 'UNKNOWN' : $seq->display_id;
  for my $f ($seq->get_SeqFeatures) 
  {
    my $ftype = $f->primary_tag;
    next unless $ftype =~ m/^($incre)$/;
    $used += $f->length;
    for my $i ($f->start .. $f->end) {
      if (exists $snp{$seqid}{$i}) {
        $count{$ftype}++;
	delete $snp{$seqid}{$i}; # dont double count!
      }
    }
  }
}

print STDERR "Calculating inter-feature SNPs.\n";
$count{'intergenic'} = $num - sum(values %count);
$count{'TOTAL'} = $num;
print Dumper(\%count, "$include used $used/$totbp bp");

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"n|noid!",  VAR=>\$noid, DEFAULT=>0, DESC=>"Process without seq IDs"},
    {OPT=>"i|include=s",  VAR=>\$include, DEFAULT=>'CDS,tRNA,rRNA', DESC=>"Include these features only, comma sep."},
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
  print "Usage: $0 [options] <reference.gbk> <nesoni/report.txt\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
