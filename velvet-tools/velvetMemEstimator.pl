#!/usr/bin/perl

$ReadSize = $ARGV[0];
$GenomeSize = $ARGV[1];
$NumReads = $ARGV[2];
$K = $ARGV[3];

$usage = "velvetMemEstimator.pl <read size> <genome size> <number of reads> <k value>\n";
$usage .= "Where:\n\t<read size> = number of base pairs in reads (20 - 100)\n";
$usage .= "\t<genome size> = estimate of number megabases in genome (2 - 10)\n";
$usage .= "\t<number of reads> = Number of reads in millions (5 - 100)\n";
$usage .= "\t<k value> = the k value used in the velveth run\n";

unless($ReadSize && $GenomeSize && $NumReads && $K){ die $usage; }


$velvetgmem = -109635 + 18977*$ReadSize + 86326*$GenomeSize + 233353*$NumReads - 51092*$K;

print "Memory estimate for velvetg: $velvetgmem bytes +/- 500Mb.\n";

$vgmGB = ($velvetgmem / 1024)/1024;

print " or $vgmGB Gb +/- 0.5Gb\n";
