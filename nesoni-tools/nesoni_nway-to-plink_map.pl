#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my(@Options, $SEP, $verbose, $prefix, $phenofile);
setOptions();

$prefix or die "Please set --prefix";

my $N;
my @S;
my $offset;
my %cds;
my %seq_of;
my %prod_of;
my %pheno;

if (-r $phenofile) {
  print STDERR "Reading: $phenofile\n";
  open PHENO, $phenofile;
  while (<PHENO>) {
    chomp;
    my @x = split ' ';
    $pheno{$x[0]} = $x[1];  
  }
  printf STDERR "Loaded %s phenotypes.\n", scalar(keys %pheno);
}

my $nwayfile = $ARGV[0];

print STDERR "Loading: $nwayfile\n";
open NWAY, $nwayfile;

while (my $line = <NWAY>) {
  chomp $line;
  my @x = split m/\t/, $line;
  if ($x[0] eq 'Reference') {
    $N = (@x - 3) / 3;
    @S = @x[ 3 .. 3+$N-1 ];
    print STDERR "Guessed $N strains: @S\n";
    $offset = 3 + 2*$N;
  }
  elsif ($line =~ m/\bCDS\b/) {
    for my $i ($offset .. $offset+$N-1) {
      my $strain = $S[$i-$offset];
      if (defined $x[$i] and $x[$i] =~ m/CDS\s+(\w+=>\w+)\s+(\S+) base \d+ codon \d+ ([^,]+)/) {
        my($diff,$gene,$product) = ($1,$2,$3);
        print STDERR "Found: $strain $gene $diff\n" if $verbose;
        $cds{$gene}->[$i-$offset]++;
        $seq_of{$gene} = $x[0];
        $prod_of{$gene} = $product;
      }
    }
  }
}
#print Dumper(\%cds);

# .map file
print STDERR "Writing: $prefix.map\n";
open MAP, ">$prefix.map";
for my $gene (sort keys %cds) {
  my $desc = $gene.'_'.$prod_of{$gene};
  $desc =~ s/\s+/_/g;
  print MAP join($SEP, $seq_of{$gene}, $desc, 0, $gene),"\n";
}

# .ped file
print STDERR "Writing: $prefix.ped\n";
open PED, ">$prefix.ped";
for my $i (0 .. $#S) {
  print PED join($SEP,
    $phenofile || $nwayfile,	# "family"
    $S[$i],     # strain
    0,0,0,      # paternal,maternal,sex
    $pheno{$S[$i]} || 'Phenotype',
  );
  for my $gene (sort keys %cds) {
#    my $changed = defined $cds{$gene}[$i] ? 1 : 0;
    my $changed = defined $cds{$gene}[$i] ? 'A' : 'T';   # Sarah email July 2014
    print PED "$SEP$changed$SEP$changed";
  }  
  print PED "\n";
}

print STDERR "Done.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"sep=s",  VAR=>\$SEP, DEFAULT=>' ', DESC=>"Output separator"},
    {OPT=>"prefix=s",  VAR=>\$prefix, DEFAULT=>'', DESC=>"Output file prefix (for .map .ped)"},
    {OPT=>"phenotypes=s",  VAR=>\$phenofile, DEFAULT=>'', DESC=>"Optional phenotype, rows of: strain pheno"},
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
  print "Usage: $0 [options] nway.table  # will output .map and .ped file\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
