#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Fatal;
use File::Spec;

my(@Options, $verbose, $source, $dir, $aligner, $threads, $ref, 
             $prefix, $draft, $closed, $binomial);
setOptions();

#----------------------------------------------------------------------

my @allid;

for my $exe (qw(date nucmer promer show-snps)) {
  require_exe($exe) ? msg("Found '$exe' - ok") : err("Can't find '$exe' in \$PATH");
}

($aligner eq 'nucmer' or $aligner eq 'promer') or err("Only nucmer/promer supported for --tool");

-r $ref or err("Please supply a readable reference .fasta with --ref");
$ref = File::Spec->rel2abs($ref); # if not File::Spec->file_name_is_absolute($ref);
msg("Reference: $ref"); 

$dir or err("Please supply a results folder with --dir");
-d $dir and err("Folder '$dir' already exists, sorry.");
msg("Making working folder: $dir");
mkdir $dir or err("Unable to make folder '$dir', sorry.");

my @denovo = map { File::Spec->rel2abs($_) } @ARGV;

chdir $dir or err("Couldn't change into '$dir', bummer.");

open MF, '>Makefile';
select MF;

print "all: ref realall\n\n";
print "ref: $ref\n\tcp -vf $ref reference.fa\n\n";
$ref = 'reference.fa';

# NC_012759.fna.gz
if ($closed) {
  my $dir = get_dirs('Bacteria', $binomial);
  for my $id (sort keys %{$dir}) {
    print "$id: \n";
    print "\tmkdir -v -p $id\n";
    print "\tzcat -f $dir->{$id}/*.fna.gz > $id/contigs.fna\n";
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
    process_it($id);
  }
}

# De novo contig sets
msg("Scanning: user supplied contigs");
for my $fn (@denovo) {
  next unless -r $fn;
  my $id = id_from_path($fn).".dir";
  print "$id:\n";
  print "\tmkdir -v -p $id\n";
  print "\tcp -v $fn $id/contigs.fna\n";
  process_it($id);
}

if (@allid <= 0) {
  err("No closed, draft or denovo genomes were found! Did you enable -c or -f ?");
}

print "realall: @allid\n";
print "\ttouch make.done\n";
#print "\trm -f genomes.txt\n";
#for my $g (@allid) {
#  print "\techo \"$g\" >> genomes.txt\n"
#}

print "clean:\n";
print "\trm -frv @allid\n";

select STDOUT;
run_cmd("nice make -j $threads");

#collate all the SNPs!
my %snp;
for my $id (@allid) {
  msg("Reading in SNPs for $id");
  open SNP, "$id/out.snps";
  while (<SNP>) {
    next unless m/^(\d+)\t([AGTC-])\t([AGTC-])\t/;
    $snp{$1}{$id} = $3;
#    msg("$id,$1,$2,$3");
  }
  close SNP;
}
msg("Found ".scalar(keys %snp)." possible SNP sites");
print Dumper(\%snp); 

my $N = @allid;
my %bases;
for my $pos (sort keys %snp) {
  my $site = $snp{$pos};
  if (scalar keys %{$site} == $N) {
    my $s = join '', values %{$site};
    msg("$s");
#    next if is_homopolymer($s);
    for my $g (sort keys %{$site}) {
      $bases{$g} .= $site->{$g};
    }
  }
}

open ALN, ">snps.faln";
for my $g (keys %bases) {
  print ALN ">$g\n",$bases{$g},"\n";
}
close ALN;
msg("Alignment written: $dir/snps.faln");

#----------------------------------------------------------------------

sub is_homopolymer {
  my($s) = @_;
  my $needle = substr $s, 0, 1;
  my $haystack = substr $s, 1;
  msg("checking $needle in $haystack");
  return 0 if $haystack =~ m/[^$needle]/;
  return 1;
}

#----------------------------------------------------------------------

sub process_it {
  my($id) = @_;
  push @allid, $id;
  msg("Got: $id");
  print "\t$aligner -p $id/out $ref $id/contigs.fna\n";
  print "\tshow-snps -C -r -T $id/out.delta > $id/out.snps\n\n";
  print "\tshow-coords -r -T $id/out.delta > $id/out.coords\n\n";
}

#----------------------------------------------------------------------

sub get_dirs {
  my($class, $binomial) = @_;
  my %d;
  msg("Scanning: $source/$class/{$binomial}");
  for my $src (qx(find $source/$class -type d)) {
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
#-------------------------------------------------------------------------

sub id_from_path {
  my($s) = @_;
  my @s = split m{/}, $s;
  $s = $s[-1];
  return $s;
}

#-------------------------------------------------------------------------

sub num_cpu {
  if ($^O =~ m/linux/i) {
    my($num) = qx(grep -c ^processor /proc/cpuinfo);
    chomp $num;
    return $num if $num =~ m/^\d+/;
  }
  return 1;
}

#-------------------------------------------------------------------------

sub run_cmd {
  msg("Running: @_");
  if (system(@_)) {
    msg("ERROR $? : $!");
    exit $?;
  }
}

#-------------------------------------------------------------------------

sub require_exe {
  my($bin) = shift;
  for my $dir (File::Spec->path) {
    my $exe = File::Spec->catfile($dir, $bin);
    return $exe if -x $exe; 
  }
  return;
}
#----------------------------------------------------------------------

sub msg {
  my($time) = qx(date +"%F %T");
  chomp $time;
  print STDERR "[$time] @_\n";
}

#----------------------------------------------------------------------

sub err {
  print STDERR "ERROR: @_\n";
  exit -1;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"s|source=s",  VAR=>\$source, DEFAULT=>'/bio/db/ncbigenomes', DESC=>"Source mirror folder"},
    {OPT=>"b|binomial=s",  VAR=>\$binomial, DEFAULT=>'', DESC=>"Pattern of which genomes to include from --source"},
    {OPT=>"d|dir=s",  VAR=>\$dir, DEFAULT=>'', DESC=>"Prefix for strain-named folders"},
    {OPT=>"c|closed!",  VAR=>\$closed, DEFAULT=>0, DESC=>"Include closed genomes"},
    {OPT=>"f|draft!",  VAR=>\$draft, DEFAULT=>0, DESC=>"Include draft genomes"},
    {OPT=>"p|prefix=s",  VAR=>\$prefix, DEFAULT=>'', DESC=>"Prefix for strain-named folders"},
    {OPT=>"a|aligner=s",  VAR=>\$aligner, DEFAULT=>'nucmer', DESC=>"Aligner: numcer or promer?"},
    {OPT=>"t|threads=i",  VAR=>\$threads, DEFAULT=>num_cpu(), DESC=>"Number of CPUs to use"},
    {OPT=>"r|ref=s",  VAR=>\$ref, DEFAULT=>'', DESC=>"Reference genome"},
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
