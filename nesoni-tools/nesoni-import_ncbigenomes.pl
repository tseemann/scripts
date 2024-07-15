#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
#use Bio::SeqIO;
#use Bio::SearchIO;

my(@Options, $verbose, $folder, $prefix, $draft, $closed, $sra, $all);
setOptions();

#----------------------------------------------------------------------

# need species etc
my $binomial = shift @ARGV or usage();

# turn --all on!
$draft=$sra=$closed=1 if $all;

#----------------------------------------------------------------------

# NC_012759.fna.gz
if ($closed) {
  my $dir = get_dirs('Bacteria', $binomial);
  print STDERR Dumper($dir) if $verbose;
  for my $id (sort keys %{$dir}) {
    print "\n# $id\n";
    print "mkdir -v -p $id\n";
    print "zcat -f $dir->{$id}/*.fna.gz > $id/contigs.fna\n";
    print "nesoni shred: --size 250 --stride 5 $id/contigs.fna > $id/s_single.fa\n";
  }
}

# NZ_AAJW00000000.contig.fna.tgz
if ($draft) {
  my $dir = get_dirs('Bacteria_DRAFT', $binomial);
  print STDERR Dumper($dir) if $verbose;
  for my $id (sort keys %{$dir}) {
    print "\n# $id\n";
    print "mkdir -v -p $id\n";
    print "tar --to-stdout -zxf $dir->{$id}/*.contig.fna.tgz > $id/contigs.fna\n";
    print "nesoni shred: --size 250 --stride 5 $id/contigs.fna > $id/s_single.fa\n";
  }
}

# Not implemented
if ($sra) {
  print STDERR "NOTICE: option --sra is not implemented yet.\n";
}

#----------------------------------------------------------------------

sub get_dirs {
  my($class, $binomial) = @_;
  my %d;
  for my $src (qx(find $folder/$class -type d)) {
    next unless $src =~ m/$binomial/i;
    chomp $src;
#    print STDERR "# $src\n";
    my @path = split m{/}, $src;
    my $label = $path[-1];
    $label =~ s/_+uid\d+$//;
#    print STDERR "# $label\n";
    $label =~ s/^[A-Z]+_[A-Z]+_+//i;
    $label =~ s/_+/_/g;
    next unless $label;
    $label = "$prefix$label" if $prefix;
    print STDERR "Found: $src => $label\n" if $verbose;
    $d{$label} = $src;
  }
  return \%d;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"f|folder=s",  VAR=>\$folder, DEFAULT=>'/bio/db/ncbigenomes', DESC=>"Mirror folder"},
    {OPT=>"p|prefix=s",  VAR=>\$prefix, DEFAULT=>'', DESC=>"Prefix for strain-named folders"},
    {OPT=>"c|closed!",  VAR=>\$closed, DEFAULT=>0, DESC=>"Include closed genomes"},
    {OPT=>"d|draft!",  VAR=>\$draft, DEFAULT=>0, DESC=>"Include draft genomes"},
    {OPT=>"s|sra!",  VAR=>\$sra, DEFAULT=>0, DESC=>"Include SRA reads [not implemented]"},
    {OPT=>"a|all!",  VAR=>\$all, DEFAULT=>0, DESC=>"Include ALL sources"},
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
  print "Usage: $0 [options] <Genus | Genus_species>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
