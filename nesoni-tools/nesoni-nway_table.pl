#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use List::MoreUtils qw(pairwise each_array);
use List::Util qw(sum);

no warnings qw(once); #Name "main::a" used only once: possible typo

my(@Options, $verbose, $prefix, $anno);
setOptions();

my $refname = "$ARGV[0]/reference.fa";
print STDERR "Loading: $refname";
my $refio = Bio::SeqIO->new(-file=>$refname, -format=>'fasta');
my $ref = $refio->next_seq->seq;
my $len = length($ref);
print STDERR " - $len bp (reference genome)\n";

my %gene_at;

if (-r $anno) {
  print STDERR "Loading: $anno";
  my $gb = Bio::SeqIO->new(-file=>$anno, -format=>'genbank');
  my $ft = $gb->next_seq;
  printf STDERR " - %d bp (reference annotation)\n", $ft->length;
  die "annotation $anno does not match $refname length" if $ft->length != $len;
  for my $f ($ft->get_SeqFeatures) {
    next if $f->primary_tag eq 'source';
    my $g = ''; ($g) = $f->get_tag_values('gene') if $f->has_tag('gene'); 
    my $p = ''; ($p) = $f->get_tag_values('product') if $f->has_tag('product');
    $p ||= $f->primary_tag;
    for my $loc ($f->location->each_Location) {
      for my $i ($loc->start .. $loc->end) {
        $gene_at{$i}{GENE} = $g;
        $gene_at{$i}{PROD} = $p;
        print STDERR "\rAnnotating: $i" if $i%101==0;
      }
    }
  }
  for my $i (1 .. $len) {
    $gene_at{$i}{GENE} ||= '';
    $gene_at{$i}{PROD} ||= '';
  }
}
print STDERR "\n";
#print Dumper($gene_at{100});

my %snp;
for my $ndir (@ARGV) {
  load_snps($ndir);
}

print STDERR "Writing: $prefix.snps.csv\n";
open my $snpfh, '>', "$prefix.snps.csv";

my @gid = (sort keys %snp);
print $snpfh tsv('Position', 'Reference', @gid, '/gene', '/product');
my $snprows=0;

my %tree;
for my $i (0 .. $len-1) {
  print STDERR "\rProcessing: $i/$len bp" if not $verbose and $i%9371==0;
  my $base = substr($ref,$i,1);
#  my $SAME = $base;
  my $SAME = '.';
  my @snp = map { $snp{$_}{$i+1} || $SAME } @gid;
  if (grep { $_ ne $SAME } @snp) {
    print $snpfh tsv($i+1, $base, @snp, $gene_at{$i}{GENE}, $gene_at{$i}{PROD} );
    $tree{'Reference'} .= $base;
    my $ea = each_array(@gid, @snp);
    while ( my ($g, $s) = $ea->() ) {
      $tree{$g} .= substr $s, 0, 1; # hack to handle indels > 1bp
    }
    $snprows++;
  }
}
print STDERR "\nWrote $snprows overall SNP rows to table.\n";

print STDERR "Writing: $prefix.tree.fasta\n";
#open my $treefh, '>', "$prefix.tree.fasta";
my $treeio = Bio::SeqIO->new(-file=>">$prefix.tree.fasta", -format=>'fasta');
my $desc = sprintf "%d SNPs across %d genomes", $snprows, 1+scalar(@gid);
for my $g (sort keys %tree) {  # should be same as (@gid) 
  $treeio->write_seq(Bio::Seq->new(
    -id=>$g, 
    -seq=>$tree{$g}, 
    -desc=>$desc,
  ));
}

#print STDERR Dumper(\%snp);

#----------------------------------------------------------------------

sub tsv {
  join(",", @_)."\n";
}

#----------------------------------------------------------------------
#Sequence        Position in reference   Change type     Old     New     Evidence
#Saa_6008        10866   substitution    T       C       "C"x39

sub load_snps {
  my($dir) = @_;
  my $repfn = "$dir/report.txt";
  -r $repfn or die "ERROR: could not read '$repfn'";
  print STDERR "Loading: $repfn";
  open my $fh, '<', $repfn;
  while (<$fh>) {
    my @x = split m/\t/;
    next unless @x==6 and $x[1] =~ m/^\d+$/;
    my $change = ($x[2] =~ m/insert/) ? lc($x[4]) : uc($x[4]);
    $snp{$dir}{$x[1]} = $change;
  }
  print STDERR " - read ", scalar(keys %{$snp{$dir}}), " SNPs\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"p|prefix=s",  VAR=>\$prefix, DEFAULT=>'nway', DESC=>"Prefix for output files: .tree.fasta snps.csv "},
    {OPT=>"a|anno=s",  VAR=>\$anno, DEFAULT=>'', DESC=>"Genbank file with annotation for reference"},
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
  print "Usage: $0 [options] nesoniDir1 [nesoniDir2 ...]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
