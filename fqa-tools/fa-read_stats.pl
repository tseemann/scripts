#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw(min max);

my $num=0;
my $sum=0;
my $min=1E99;
my $max=0;
my $len=0;
my @freq;

my $file = @ARGV ? join(';', @ARGV) : "stdin";

while (<ARGV>) {
  if (m/^>/) {
    $num++;
    next if $num==1; # no stats yet
    $sum += $len;
    $min = min($min, $len);
    $max = max($max, $len);
    $freq[$len]++;
    $len = 0; # reset for new sequence
  }
  else {
    chomp;
    $len += length;
  }
}

my $tot=0;
my $median = $min;
while ($tot <= $num/2) {
  $tot += ($freq[$median] || 0);
  $median++;
}
$median--;

my $mean = int($sum / $num);
print "FILE\tNSEQ\tBASES\tMIN\tMEAN\tMEDIAN\tMAX\n";
print "$file\t$num\t$sum\t$min\t$mean\t$median\t$max\n";
