#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum min max);

for my $file (@ARGV) {
  process($file);
}

sub process {
  my $fn = shift;
  open FASTA, $fn or return;
  my $n=0;
  my $bp = -1;
  my $Ns = 0;
  my %size;
  while (<FASTA>) {
    chomp;
    if (m/>/) {
      $n++;
      if ($bp >= 0) {
        $size{$bp}++;
      }
      $bp=0;
    }
    else {
      $bp += length;
      s/[^N]//gi;
      $Ns += length;
    }
  }     
  
  my $total = total(\%size);
  my $mean = $total / $n;
  
  printf "file=%s nseq=%d total=%d Ns=%d min=%d mean=%d med=%d max=%d N50=%d\n", 
    $fn, 
    $n, 
    $total,
    $Ns,
    min(keys %size),
    $mean,
    -1,
    max(keys %size),
    n50($total, \%size),
    ;
}

sub total {
  my $h = shift;
  return sum( map { $h->{$_} * $_ } keys %$h );
}

sub n50 {
  my($total, $h) = @_;
  my $sum=0;
  $total = $total/2;
  for my $len (keys %$h) {
    $sum += $h->{$len} * $len;
    return $len if $sum >= $total;
  }
  return -1;
}

sub mean {
  my($n, $h) = @_;
  my $sum=0;
  for my $len (keys %$h) {
    $sum += $h->{$len} * $len;
  }
  return $sum/$n;
}
  
