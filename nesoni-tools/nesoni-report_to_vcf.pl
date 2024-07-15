#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use Time::Piece;
use List::Util qw(sum max);
use File::Spec;

my(@Options, $verbose);
setOptions();

for my $ndir (@ARGV) {
  print STDERR "Processing: $ndir\n";

  my %snp;
  my $report = "$ndir/report.txt";
  next unless -r $report;
  open REPORT, '<', $report;
  while (<REPORT>) {
    chomp;
    my @x = split m/\t/;
    next unless @x > 1 and $x[1] =~ m/^\d+$/;
    # sum up evidence to get overall depth
    my $depth = $x[5];
    $depth =~ s/\D+/ /g;
    my @d = split ' ', $depth;
    $depth = sum(@d);
    # estimate quality
    my $qual = int (-10 * log(1 - (max(@d) / ($depth+1)) ) / log(10) );
    # save for later
    $snp{ $x[0] }{ $x[1] } = [ $x[3], $x[4], $depth, $qual ];
  }
  close REPORT;
  print Dumper(\%snp) if $verbose;

  my %seq;
  my $ref = "$ndir/reference/reference.fa";
  next unless -r $ref;
  my $sin = Bio::SeqIO->new(-file=>$ref, -format=>'fasta');
  while (my $seq = $sin->next_seq) {
    my $id = $seq->id;
    die "Sequence ID '$id' not allowed to have ':' in it in VCF" if $id =~ m/:/;
    $seq{$id} = $seq->seq;
  }
  
  my $vcf = "$ndir/report.vcf";
  open VCF, '>', $vcf;
  select VCF;

  print "##fileformat=VCFv4.1\n";
  my $t = localtime;
  print "##fileDate=", $t->ymd(''), "\n";
  print "##source=$0\n";
  printf "##reference=file://%s\n", File::Spec->rel2abs($ref);
  print "##INFO=<ID=DP,Number=1,Type=Integer,Description='Total Depth'>\n";
  for my $id (sort keys %seq) {
    printf "##contig=<ID=$id,length=%d>\n", length($seq{$id});
  }
  print "#", join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO)),"\n";
  for my $id (sort keys %snp) {
    for my $pos (sort { $a <=> $b } keys %{$snp{$id}} ) {
      my($ref,$alt,$depth, $qual) = @{ $snp{$id}{$pos} };
      #print "$ref => $alt\n";
      if ($alt eq '-') {
        $pos--;
        my $prev = substr $seq{$id}, $pos, 1;
	$alt = $prev;
	$ref = $prev.$ref;
      }
      elsif ($ref eq '-') {
        $pos++;
        my $prev = substr $seq{$id}, $pos, 1;
	$alt = $prev.$alt;
	$ref = $prev;
      }
      print "$id\t$pos\t.\t$ref\t$alt\t$qual\tPASS\tDP=$depth\n";
    }
  }
  close VCF;
  
  #system("bgzip -c $vcf > ${vcf}.gz");
  system("bgzip -f $vcf");
  system("tabix -f -p vcf ${vcf}.gz");
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
  print "Usage: $0 [options] <nesoni_consensus_dir1> ...\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
