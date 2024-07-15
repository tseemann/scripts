#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use Algorithm::Cluster;
#use Tree::Simple;
use SVG;

my(@Options, $verbose, $vsize, $cluster);
setOptions();

# Species       Genes   Alg.-Conn.      FP929058.faa    CP004063.faa    CP006620.faa ...

my @species;
my $S;
my @cluster;
my $N;
my @matrix;

print STDERR "Reading .proteinortho\n";

while (<ARGV>) {
  chomp;
  my @x = split m/\t/;
  if ($x[0] =~ m/^# Species/) {
    $S = @x - 3;
    @species = splice @x, 3;
    die "bad S calc" if $S != @species;
  }
  else {
    push @cluster, [ @x ];
    push @matrix, [ map { $_ eq '*' ? 0 : 1 } (splice @x, 3) ];
  }
}
$N = @cluster;

#print STDERR Dumper(\@matrix);
my $tree = Algorithm::Cluster::treecluster(
  data => \@matrix,
  mask => '',
  weight => '',
  transpose => 1,
  dist => 'c',
  method => 'a',
);
$tree->scale;
#print STDERR Dumper($tree);
print STDERR "S=$S and treeLength=",$tree->length,"\n";
for (my $i=0; $i < $tree->length; $i++) {
  my $node = $tree->get($i);
  printf STDERR "%15s [%3d] %3d %3d %7.3f\n",$species[$i],-1-$i,$node->left,$node->right,$node->distance;
#  print STDERR Dumper($tree->get($i));  
}

my %t;
my %has_parent;

for (my $i=0; $i < $tree->length; $i++) {
  my $node = $tree->get($i);
  my $id = -1-$i;
  $t{$id}{L} = $node->left;
  $t{$id}{R} = $node->right;
  $has_parent{$node->left}++ if $node->left < 0;
  $has_parent{$node->right}++ if $node->right < 0;
}
my $root = 0;
for (my $i=0; $i < $tree->length; $i++) {
  if (not $has_parent{-1-$i}) {
    $root = -1-$i;
    last;
  }
}
my @order;

sub inorder {
  my $id = shift;
  # left
  my $value = $t{$id}{L};
  if ($value < 0) {
    inorder($value);
  } 
  else {
    push @order, $value;
  }
  # right
  $value = $t{$id}{R};
  if ($value < 0) {
    inorder($value);
  } 
  else {
    push @order, $value;
  }
}

inorder($root);
@order = reverse @order;

print STDERR Dumper(\%t, \%has_parent, $S, $root, \@order);

#exit;

my $sort_col = 3;
print STDERR "Sorting $N clusters by column $sort_col\n";
@cluster = reverse sort { $a->[$sort_col] cmp $b->[$sort_col] } @cluster;

print STDERR "Drawing SVG\n";
my $svg = SVG->new(width=>$N, height=>(2*$S+1)*$vsize);

my $fontsize = int(0.5*$vsize);
my $x = 0;
for my $c (@cluster) {
  for my $y (1 .. $S) {
    my $s = $order[$y-1];
    my $colour = ($c->[3+$s] eq '*') ? 'lightgray' : 'green';
    $svg->rectangle(
      x=>$x, y=>(($y-1)*$vsize+$y), width=>1, height=>$vsize,
      style=>{ fill=>$colour }
    );
    if ($x==$N-1) {
      $svg->text( 
        x=>0, y=>$y*($vsize+1) - int(0.2*$vsize),
        style=>{ 
          'fill'=>'black',
          'text-anchor'=>'start', 
          'font-size'=>"${fontsize}px", 
          'font-family'=>'sans-serif',
        } 
      )->cdata($species[$s]);
    }
  }
  $x++;
}

print $svg->xmlify;

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"vsize=i",  VAR=>\$vsize, DEFAULT=>20, DESC=>"Species height"},
    {OPT=>"cluster!",  VAR=>\$cluster, DEFAULT=>0, DESC=>"Hierarchial clustering of species"},
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
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
