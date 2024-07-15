#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
#use Bio::SeqIO;
#use Bio::SearchIO;

my(@Options, $verbose, $folder, $prefix, $draft, $closed, $sra, $all, $binomial);
setOptions();

#----------------------------------------------------------------------

# turn --all on!
$draft=$sra=$closed=1 if $all;

#----------------------------------------------------------------------

my $ref;
my @allid;

print "all: realall\n";

# NC_012759.fna.gz
if ($closed) {
  my $dir = get_dirs('Bacteria', $binomial);
  for my $id (sort keys %{$dir}) {
    print "$id:\n";
    print "\tmkdir -v -p $id\n";
#    print "\tzcat -f $dir->{$id}/*.fna.gz > $id/contigs.fna\n";
    print "\tcat $dir->{$id}/*.fna > $id/contigs.fna\n";
    $ref ||= "$id/contigs.fna";
    process_it($id);
  }
}

# NZ_AAJW00000000.contig.fna.tgz
if ($draft) {
  my $dir = get_dirs('Bacteria_DRAFT', $binomial);
  for my $id (sort keys %{$dir}) {
    print "$id:\n";
    print "\tmkdir -v -p $id\n";
    print "\trm -f $id/contigs.fna\n";
    for my $tarball ( glob("$dir->{$id}/*.contig.fna*tgz") ) {
      print "\ttar --to-stdout -zxf $tarball >> $id/contigs.fna\n";
    }
    $ref ||= "$id/contigs.fna";
    process_it($id);
  }
}

# Not implemented
if ($sra) {
  print STDERR "NOTICE: option --sra is not implemented yet.\n";
}

# De novo contig sets
for my $fn (@ARGV) {
  my @id = split m{/}, $fn;
  my $id = $id[-1];
  print "$id:\n";
  print "\tmkdir -v -p $id\n";
  print "\tcp '$fn' $id/contigs.fna\n";
  $ref ||= "$id/contigs.fna";
  process_it($id);
}

print "realall: @allid\n";
print "\trm -f genomes.txt\n";
for my $g (@allid) {
  print "\techo \"$g\" >> genomes.txt\n"
}
print "\tnesoni nway: --as nexus --reference no --require-all yes --indels no @allid > tree.nex\n";

#----------------------------------------------------------------------

sub process_it {
  my($id) = @_;
  push @allid, $id;
  print STDERR "Processing: $id\n";
  print "\tnesoni shred: --size 50 --stride 5 $id/contigs.fna > $id/s_single.fa\n";
  print "\tnice nesoni shrimp: $id --sam-unaligned no $ref reads: $id/s_single.fa shrimp-options: --threads 1 --local\n";
  print "\tnice nesoni consensus: $id --cutoff 0.1\n\n";
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
    $label =~ s/serovar.//;
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
    {OPT=>"b|binomial=s",  VAR=>\$binomial, DEFAULT=>'', DESC=>"Pattern of which genomes to include"},
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
  print "Usage: $0 [options] <contigs1.fna ...>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
