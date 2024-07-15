#!/usr/bin/env perl
use strict;
use IO::Handle;
use Data::Dumper;
use List::Util qw(max);

my(@Options, $verbose, $K, $C);
setOptions();

my $in = IO::Handle->new_from_fd( fileno(STDIN), 'r' );
my $peek = $in->getc();
print STDERR "First char is '$peek'.\n";
$in->ungetc(ord($peek));
my $nline = $peek eq '>' ? 2 : 4;
print STDERR "Each sequence uses $nline lines.\n";

my %freq;
my $count=0;

while (not $in->eof) {
  my @s = map { $in->getline } (1 .. $nline);
#  print $s[1];
  chomp $s[1];
  kmers(uc $s[1]);
  $count++;
  print STDERR "\rProcessed reads: $count" if $count % 1E4 == 0;
}
print STDERR "\n";
print STDERR "Hash has ", scalar(keys %freq), " k-mers.\n";

print Dumper(\%freq);

my %num;
map { $num{$_}++ } values %freq;
print "$_,$num{$_}\n" for (sort { $a<=>$b } keys %num);

#print Dumper(\%num);  

#----------------------------------------------------------------------

sub kmers {
  my $s = shift;
  my $L = length($s);
  for (my $i=0; $i < $L-$K+1; $i++) {
    $freq{ substr($s, $i, $K) } ++;
  }
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"k", VAR=>\$K, DEFAULT=>21, DESC=>"k-mer size"},
    {OPT=>"c", VAR=>\$C, DEFAULT=>4, DESC=>"Chunk size for data structure"},
  );

#  (@ARGV < 2) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 {file.fa | file.fq}\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
