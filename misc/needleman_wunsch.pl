#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw(max);
use Data::Dumper;

my @in = (
  [ 'A', 'T' ],
  [ 'ATGC', 'CTATAG' ],
  [ 'A', 'T' ],
  [ 'ATC', 'ACT' ],
  [ 'AAA', 'ATAA' ],
  [ 'AGGT', 'AGCT' ],
  [ 'AGGT', 'AGGCT' ],
  [ 'AGGT', 'AG' ],
  [ 'AGGT', 'AGGGT' ],
  [ 'AGGT', 'AGGGGGT' ],
);
  
for my $in (@in) {
  print "\nIN:\n$in->[0]\n$in->[1]\n";
  my @out = needleman_wunsch(@$in);
  print "\nOUT:\n$out[0]\n$out[1]\n";
}

sub needleman_wunsch {
  my($x, $y, $gap) = @_;

  $gap ||= -2;

  sub S {
    my($c1,$c2) = @_;
    my $score = ($c1 eq $c2 ? +1 : -1);
    return $score;
  }

  my $Lx = length $x;
  my $Ly = length $y;
  my @x = ('N', split m//, $x);
  my @y = ('N', split m//, $y);
  my @F;

  $F[$_][0] = $_*$gap for (0 .. $Ly); # row: y / j
  $F[0][$_] = $_*$gap for (0 .. $Lx); # col: x / i

  for my $j (1 .. $Ly) {
    for my $i (1 .. $Lx) {
      $F[$j][$i] = max(
        $F[$j-1][$i-1] + S($x[$i], $y[$j]),
        $F[$j-1][$i  ] + $gap,
        $F[$j  ][$i-1] + $gap,
      );
    }
  }
  
#  print Dumper(\@x, \@y, \@F);
  
  my $i = $Lx;
  my $j = $Ly;
  
  my $ax = '';
  my $ay = '';
  
  while ($i > 0 or $j > 0) {
#    print "i=$i j=$j ax=$ax ay=$ay x=$x y=$y\n";
#    print "Fji: .=$F[$j][$i] \\=$F[$j-1][$i-1] -=$F[$j][$i-1] |=$F[$j-1][$i]\n";
    if ($j > 0 and $F[$j][$i] == $F[$j-1][$i] + $gap) {
      $ax =    '-' . $ax;
      $ay = $y[$j] . $ay;
      $j--;
    }
    elsif ($i > 0 and $F[$j][$i] == $F[$j][$i-1] + $gap) {
      $ax = $x[$i] . $ax;
      $ay =    '-' . $ay;
      $i--;
    }
    elsif ($i > 0 and $j > 0 and $F[$j][$i] == $F[$j-1][$i-1] + S($x[$i], $y[$j]) ) {
      $ax = $x[$i] . $ax;
      $ay = $y[$j] . $ay;
      $i--;
      $j--;
    }
    else {
      die "should not be here i=$j j=$j ax=$ax ay=$ay x=@x y=@y";
    }
  }

  return ($ax, $ay);
}


