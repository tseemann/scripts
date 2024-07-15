#!/usr/bin/env perl
use warnings;
use strict;
use List::Util;
use List::MoreUtils qw(all any);
use Data::Dumper;

my(@Options, $verbose, $po_fn, $desc_fn, $inc_fn, $SEP, $ignore_fn);
setOptions();

# --ortho  is mandatory
my $po_fh = get_fh($po_fn, '--ortho');

# --include is mandatory
my @inc;
@inc = load_array( get_fh($inc_fn, '--include') ); 
printf STDERR "Using %d positive phenotype taxa: @inc\n", 0+@inc;
my %inc = (map { $_ => 1 } @inc);
print STDERR Dumper(\%inc) if $verbose;

# --ignore is optional
my @ignore;
if ($ignore_fn) {
  @ignore = load_array( get_fh($ignore_fn, '--ignore') );
  printf STDERR "Ignoring %d taxa: @ignore\n", 0+@ignore;
}
my %ignore = (map { $_ => 1 } @ignore);
print STDERR Dumper(\%ignore) if $verbose;

my $desc;
if ($desc_fn) {
  $desc = load_hash( get_fh($desc_fn, '--descs') );
  printf STDERR "Loaded %d gene descriptions.\n", scalar(keys %$desc);
}

my @exc;
my @taxa;
my %index_of;
my $N;
my $found=0;
while (<$po_fh>) {
  chomp;
  s/,/ /g; # CRITICAL LINE 3
  my @x = split m/\t/;
  if ($x[1] eq 'Genes') {
    # 1st header line
    @taxa = @x[3 .. $#x];
    $N = scalar @taxa;
    print STDERR "Found $N taxa: @taxa\n";
    %index_of = (map { $taxa[$_] => $_ } 0 .. $#taxa);
    # input checking
    my %taxon = (map { $_ => 1 } @taxa);
    foreach (@inc) {
      die "unknown taxon '$_' in $inc_fn" if not exists $taxon{$_};
    }  
    foreach (@ignore) {
      die "unknown taxon '$_' in $ignore_fn" if not exists $taxon{$_};
    }  
    @exc = grep { !$inc{$_} && !$ignore{$_} } @taxa; # CRITICAL LINE #1
    printf STDERR "Ignoring %d taxa: @ignore\n", scalar @ignore;
    printf STDERR "Using %d +ve taxa: @inc\n", scalar @inc;
    printf STDERR "Using %d -ve taxa: @exc\n", scalar @exc;
    # output header
    print join($SEP,'COUNTER', 'DESCRIPTION', @inc),"\n";
  }
  else {
    # cluster (gene ortholog) line
    my %have;
    my @gene = @x[3 .. $#x];
    for my $i (0 .. $N-1) {
      #print "$i $taxa[$i] $gene[$i]\n";
      if ($gene[$i] ne '*') {
        $have{ $taxa[$i] } = $gene[$i];
      }
    }
    print STDERR Dumper(\%have, \@inc) if $verbose;
    next if any { $have{$_} } @exc; # CRITICAL LINE #2
    my @col =  map { $gene[$index_of{$_}] } @inc;
    my @notempty = grep { $_ ne '*' and $_ !~ m/,/ } @col;
    my $label = $desc->{ $notempty[0] } || 'no description, use --descs';
    $label =~ s/,/ /g;
    $found++;
    print join($SEP, $found, $label, @col),"\n";
  }  
}
print STDERR "Found $found matches.\n";

#----------------------------------------------------------------------

sub get_fh {
  my($fname, $desc) = @_;
  $fname or die "Please specify file for '$desc'.";
  -r $fname or die "File '$fname' is not readable.";
  open my $fh, '<', $fname or die "Could not open file '$fname'";
  print STDERR "Opened $desc file: $fname\n";
  return $fh;
}

#----------------------------------------------------------------------

sub load_array {
  my($fh) = @_;
  my @x = <$fh>;
  chomp @x;
  @x = grep { length($_) > 0 } @x;
  return @x;
}

#----------------------------------------------------------------------

sub load_hash {
  my($fh, $sep) = @_;
  $sep ||= qr/\t/;
  my %res;
  while (<$fh>) {
    chomp;
    my @x = split $sep;
    $res{ $x[0] } = $x[1];
  }
  return \%res;
}

#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"separator=s",  VAR=>\$SEP, DEFAULT=>"\t", DESC=>"Output separator (default=<TAB>)"},
    {OPT=>"ortho=s",  VAR=>\$po_fn, DEFAULT=>'', DESC=>"Input .proteinortho file"},
    {OPT=>"descs=s",  VAR=>\$desc_fn, DEFAULT=>'', DESC=>"Input .descriptions file (optional)"},
    {OPT=>"include=s",  VAR=>\$inc_fn, DEFAULT=>'', DESC=>"File with taxa with phenotype"},
    {OPT=>"ignore=s",  VAR=>\$ignore_fn, DEFAULT=>'', DESC=>"File with taxa to ignore (optional)"},
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
  print "Usage: $0 [options] --ortho PAN.proteinortho --desc PAN.descriptions --include INCLUDE.TXT > OUTPUT.csv\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__

Margaret and I want to look at CDS unique to this cluster of isolates (ST796):

Ef_aus1154
Ef_aus1139
Ef_aesops12018
Ef_aus1095
Ef_aus1173
with the option of including Ef_aesops13004 (this isolate has somewhat lower isoprop tolerance than its other ST796 colleagues).

Our pan genome analysis is here:
http://dna.med.monash.edu.au/~tstinear/fripan/pan.html?Ef_pan-3

This was prepared with the following command:
nice proteinortho5.pl -e=1e-09 -project=Ef_pan-3 -identity=80 -cov=80 -clean *faa &

using data found here:
/home/tstinear/Projects.dir/Ef.dir/Ef.map/Ef_12018.map/BratNextGen/faa-files

Thanks

Tim
