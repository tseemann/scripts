#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::SearchIO;

my(@Options, $verbose, $clen);
setOptions();

my $dir = shift @ARGV;
die "bad nesoni folder '$dir'" unless -d $dir;

my $rep = "$dir/report.txt";
die "SNP report '$rep' not found!" unless -r $rep;

my $ref = "$dir/reference.gbk";
$ref = shift(@ARGV) if not -r $ref;
die "reference '$ref' not found!" unless -r $ref;

my %snp;
print STDERR "Loading: $rep\n";
open SNP, '<', $rep;
my $hdr = <SNP>;
chomp $hdr;
while (<SNP>) {
  chomp;
  my($sid, $pos) = split m/\t/;
  next unless $pos =~ m/^\d+$/;
  $snp{$sid}{$pos} = $_;
}

print STDERR "Loading: $ref\n";
my $gbk = Bio::SeqIO->new(-file=>$ref, -format=>'genbank');
print "$hdr\tContext\tLocus\tGene\tProduct\n";
while (my $seq = $gbk->next_seq) {
  my $sid = $seq->display_id;
#  print STDERR "\t$sid\n";
  for my $f ($seq->get_SeqFeatures) {
    next unless $f->primary_tag eq 'CDS';
    my %col;
    for my $tag ('locus_tag', 'gene', 'product') {
      $col{$tag} = '';
      if ( $f->has_tag($tag) ) {
        ($col{$tag}) = $f->get_tag_values($tag);
      }
    }
    for my $pos ($f->start .. $f->end) {
      my $tsv = $snp{$sid}{$pos} or next;
      my $context = lc $seq->subseq($pos-$clen, $pos+$clen);
      substr($context, $clen, 1) = uc substr($context, $clen, 1);
      print "$tsv\t$context\t$col{locus_tag}\t$col{gene}\t$col{product}\n";
    } 
  }
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"c|context=i",  VAR=>\$clen, DEFAULT=>5, DESC=>"SNP context size"},
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
  print "Usage: $0 [options] <nesoni_dir> [ref.gb]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
