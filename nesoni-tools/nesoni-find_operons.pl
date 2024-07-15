#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::Sam;

my(@Options, $verbose, $mincov, $testmode);
setOptions();

my $gbfile = shift @ARGV;
-r $gbfile or die "can not read genbank file: $gbfile";
my $nesdir = shift @ARGV;
-d $nesdir or die "can not see nesoni dir: $nesdir";
#my $bamfile = "$nesdir/alignments_filtered_sorted.bam";
my $bamfile = "$nesdir/alignments.bam";
-r $bamfile or die "can not find BAM file: $bamfile";

my $SEQID = 'unknown';

my %prod_of;
my %feat_of;

my @map;
print STDERR "Loading GBK: $gbfile\n";
my $C=0;
my $gbk = Bio::SeqIO->new(-file=>$gbfile, -format=>'genbank');
while (my $seq = $gbk->next_seq) {
  $SEQID = $seq->display_id;
  print STDERR "Found: $SEQID\n";
  for my $f ($seq->get_SeqFeatures) {
    next if $f->primary_tag ne 'CDS';
    $prod_of{ tagval($f) } = tagval($f, 'product');
    $feat_of{ tagval($f) } = $f;
    my $strand = $f->strand > 0 ? 0 : 1;
    for my $i ($f->start .. $f->end) {
#      push @{$fwd[$i]}, $f;
      $map[$strand][$i] = $f;
    }
    $C++;
  }
  last; # only process FIRST seq
}
#print Dumper(\@fwd);
#print STDERR Dumper(\%prod_of); exit;
print STDERR "Number of CDS: $C\n"; 

print STDERR "Processing BAM: $bamfile\n";
my %span;
open my $sam, "samtools view $bamfile |";
while (<$sam>) {
  print STDERR "\rRead: $." if $. % 19357 == 0;
  last if $testmode and $. > 1E6 ;
  my @x = split m/\t/;
#  next unless $x[5] =~ m/^(\d+)M/; # note no '$' at end
#  my $len = $1;
  my $len = length $x[9];
  my $p1 = $x[3];
  my $p2 = $p1 + $len - 1;
#  my $strand = ($x[1] & 0x0010) ? -1 : +1;
  my $strand = ($x[1] & 0x0010) ? 1 : 0;
#  print STDERR $_;
  my $L = $map[$strand][$p1] or next;
  my $R = $map[$strand][$p2] or next;
  if ($L ne $R) {
#    print STDERR "Span: $p1..$p2\t$strand\t${len}bp\n";
#    print "\t", desc($L), "\n";
#    print "\t", desc($R), "\n";
    my $key = tagval($L).'~~~'.tagval($R);
    $span{$strand}{$key}++;
  }
}
close $sam;
print STDERR "\n";

my @list;
my $tablefn = "$nesdir/operon-pairs.csv";
print STDERR "Writing: $tablefn\n";
open my $TFH, '>', $tablefn or die "can't write to: $tablefn";
print $TFH join("\t", qw(GENE1 GENE2 STRAND NUMREADS PROD1 PROD2)),"\n";
for my $strand (0, 1) {
  for my $pair (sort keys %{$span{$strand}}) {
    my $cov = $span{$strand}{$pair};
    next unless $cov >= $mincov;
    my @pair = split m/~~~/, $pair;
    my @row = (
      @pair, 
      ($strand ? '-' : '+'),
      $cov,
      $prod_of{$pair[0]},
      $prod_of{$pair[1]}
    );
    printf $TFH join("\t",@row)."\n";
    push @list, \@row;
  }
}
close $TFH;

printf STDERR "Operon pairs: %d\n", scalar @list;
@list = chain_pairs(@list);
printf STDERR "Operon chains: %d\n", scalar @list;

my $gfffn = "$nesdir/operon-pairs.gff";
print STDERR "Writing: $gfffn\n";
open my $GFH, '>', $gfffn or die "can't write to: $gfffn";
print $GFH "##gff-version 3\n";
my $count=0;
for my $row (@list) {
  $count++;
  print $GFH join("\t",
    $SEQID,
    '.',
    'operon',
    $feat_of{$row->[0]}->start,
    $feat_of{$row->[1]}->end,
    '.',
    $row->[2],
    '.',
    "ID=txop_$count;colour=200 200 50;Note=$row->[4]"
  ), "\n";
}

#----------------------------------------------------------------------

sub chain_pairs {
  my(@list) = @_;
  my(@new) = shift @list;
  while (my $next = shift @list) {
    if ($new[$#new]->[1] ne $next->[0]) {
      push @new, $next;
    }
    else {
      $new[$#new]->[1] = $next->[1];
      print STDERR "Chaining: ",$next->[0],"\n";
    }
  }  
  return @new;
}


#----------------------------------------------------------------------

sub tagval {
  my($f, $tag) = @_;
  $tag ||= 'locus_tag';
  my $val = "no-$tag";
  if ($f->has_tag($tag)) {
    ($val) = $f->get_tag_values($tag);
  }
  return $val;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"m|mincov=i",  VAR=>\$mincov, DEFAULT=>1, DESC=>"Minimum read coverage"},
    {OPT=>"t|test!",  VAR=>\$testmode, DEFAULT=>0, DESC=>"Test mode, only examine 1E6 reads"},
  );

  (@ARGV > 1) or usage();

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] <reference.gbk> <nesoni_dir>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
