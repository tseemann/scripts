#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use List::Util qw(sum);

my(@Options, $verbose, $noid);
setOptions();

my $fisher_fn = shift @ARGV;
my $ref_fn = shift @ARGV;

print STDERR "Loading fisher report: $fisher_fn\n";
open FISHER, '<', $fisher_fn;
my %snp;
my $snps=0;
while (<FISHER>) {
  chomp;
  my @x = split m/\t/;
  #NC_012470.1     121533  substitution    C       Tx198 Cx49 Ax1  Cx193   1.31787e-78
  next unless @x == 7 and $x[1] =~ m/^(\d+)$/;
  $x[0] = 'UNKNOWN' if $noid;
  $x[0] =~ s/\.\d+$//; #hack
  $snp{$x[0]}{$x[1]} = \@x;
  $snps++;
}
print STDERR "Read $snps SNPs.\n";
die "No SNPs to test" if $snps <= 0;
print Dumper(\%snp) if $verbose;

print STDERR "Loading reference: $ref_fn\n";
my $gbk = Bio::SeqIO->new(-file=>$ref_fn, -format=>'genbank');
while (my $seq = $gbk->next_seq) {
  my $id = $noid ? 'UNKNOWN' : $seq->accession_number;
  print STDERR "Checking: $id\n";
  next unless exists $snp{$id};
  for my $f ($seq->get_SeqFeatures) {
    next unless $f->primary_tag eq 'CDS';
    my @anno;
    for my $t (qw(locus_tag gene product)) {
      push @anno, $f->has_tag($t) ? ($f->get_tag_values($t))[0] : '';      
    }
    for my $i ($f->start .. $f->end) {
      my $row = $snp{$id}{$i};
      next unless $row;
      print join("\t", @$row, @anno),"\n";
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
    {OPT=>"n|noid!",  VAR=>\$noid, DEFAULT=>0, DESC=>"Process without seq IDs"},
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
  print "Usage: $0 [options] <fisher.csv> <reference.gbk>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
