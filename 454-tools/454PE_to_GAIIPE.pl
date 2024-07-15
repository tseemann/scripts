#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use List::Util qw(min max);

my(@Options, $verbose, $minlen, $maxlen, $maxhomo, $overlap);
setOptions();
$maxhomo++;

my $counter=0;
my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');

while (my $left = $in->next_seq) {
  my $right = $in->next_seq;

  my @L = homo_split($left->seq);
  @L = grep { length($_) >= $minlen } @L;
  print STDERR "Lshort = @L\n" if $verbose;
  next unless @L > 0;

  my @R = homo_split($right->seq);
  @R = grep { length($_) >= $minlen } @R;
  print STDERR "Rshort = @R\n" if $verbose;
  next unless @R > 0;

  if (@L > @R ) {
    @R = expand_out(\@R, scalar(@L));
  }
  elsif (@R > @L) {
    @L = expand_out(\@L, scalar(@R));
  }  
  die "R and L not same size arrays" unless @R == @L;
  
  my $id = $left->display_id;  
  for my $i (0 .. $#R) {
    next unless length($L[$i]) >= $minlen and length($R[$i]) >= $minlen;
    print ">$id.$i/1\n$L[$i]\n>$id.$i/2\n$R[$i]\n";
  }

#  my $subcounter=0;
#  $counter++;
#  for my $l (@L) {
#    for my $r (@R) {
#      my $id = sprintf "%s_%d_%d", $left->display_id, $counter, ++$subcounter;
#      printf ">%s/1\n%s\n>%s/2\n%s\n", $id, $l, $id, $r;
#    }
#  }

  print STDERR "\n" if $verbose;
} 


#----------------------------------------------------------------------

sub expand_out {
  my($in, $size) = @_;
  # put in array
  my @in = @$in;
  # while we need to keep adding more sequences;
  while (@in != $size) {
    # sort so longest at start
    @in = sort { length($b) <=> length($a) } @in;
    # chop the longest into overlapping halves
    my $long = shift @in;
    my $x = substr $long, 0, (length($long)-$overlap)/2 + $overlap;
    my $y = substr $long, (length($long)-$overlap)/2;
    print STDERR "EXPAND\n$long\n$x\n\t\t$y\n\n" if $verbose;
    push @in, $x, $y;
  }
  return @in; 
}

#----------------------------------------------------------------------
# FIXME - this could be done with a simple split w/ capture 
# eg. 
#  @s = split m/(.{4,)/;  
# @s = map { $s[$_] = substr($s[$_],0,4) } (1,3,5..)
# maybe?

sub homo_split {
  my($s) = @_;
  $s = uc $s;
  print STDERR "before = $s\n" if $verbose;
  $s =~ s/A{$maxhomo,}/'a'x($maxhomo-1).'N'.'a'x($maxhomo-1)/ge;
  $s =~ s/T{$maxhomo,}/'t'x($maxhomo-1).'N'.'a'x($maxhomo-1)/ge;
  $s =~ s/G{$maxhomo,}/'g'x($maxhomo-1).'N'.'a'x($maxhomo-1)/ge;
  $s =~ s/C{$maxhomo,}/'c'x($maxhomo-1).'N'.'a'x($maxhomo-1)/ge;
#  $s =~ s/(.){$maxhomo,}/"$1"x($maxhomo-1).'N'."$1"x($maxhomo-1)/ge;
  print STDERR "middle = $s\n" if $verbose;
  my @s = split m/N+/, $s;
  print STDERR "finish = @s\n" if $verbose;
  return @s;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"m|minlen=i",  VAR=>\$minlen, DEFAULT=>29, DESC=>"Minimum read length"},
#    {OPT=>"l|maxlen=i",  VAR=>\$maxlen, DEFAULT=>100, DESC=>"Maximum/target read length"},
    {OPT=>"o|overlap=i",  VAR=>\$overlap, DEFAULT=>31, DESC=>"Overlap magic value"},
    {OPT=>"p|maxhomo=i",  VAR=>\$maxhomo, DEFAULT=>3, DESC=>"Maximum homopolymer tolerated"},
  );

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] 454_PE.fasta > GAII_PE.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
