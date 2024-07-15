#!/usr/bin/env perl
use strict;
use threads;
#use Bio::SeqIO;
use File::Temp qw(tempdir);
use Data::Dumper;

my(@Options, $verbose, $identity, @set);
setOptions();

foreach (1, 2) {
  $set[$_] or err("need protein fasta file for set $_ via option --1");
  -r $set[$_] or err("can't read file: $set[$_]");
}
#$set2 or err("need protein fasta file for set 2 via option --2");
#-r $set2 or err("can't read file: $set2");

# http://weizhong-lab.ucsd.edu/cd-hit/wiki/doku.php?id=cd-hit_user_guide
#-n 5 for thresholds 0.7 ~ 1.0
#-n 4 for thresholds 0.6 ~ 0.7
#-n 3 for thresholds 0.5 ~ 0.6
#-n 2 for thresholds 0.4 ~ 0.5 
($identity > 0.4 && $identity <= 1.0) or err("identity must be between 0.4 and 1.0");
my $wordsize = 5;
$wordsize = 4 if $identity < 0.7;
$wordsize = 3 if $identity < 0.6;
$wordsize = 2 if $identity < 0.5;

my $dir = tempdir(CLEANUP=>1);

msg("First set: $set[1]");
msg("Second set: $set[2]");
msg("Identity: $identity");
msg("Wordsize: $wordsize");
msg("Temp dir: $dir");

my %desc;
for my $s (1,2) {
  msg("Copying and re-labelling set $s");
  system("sed 's/^>/>$s~/' < $set[$s] >> $dir/proteins")==0 or err("Failed to run: sed");
  msg("Saving fasta descriptions from set $s for later");
  open FAA, '<', $set[$s];
  while (<FAA>) {
    next unless m/^>(\S+)\s+(.*)$/;
    chomp;
    $desc{$s}{$1} = $2;
  }
  close FAA;
}
#system("grep '^>' $dir/proteins"); exit;
#print Dumper(\%desc);  exit;

msg("Cluster proteins with cd-hit");
my $cmd = "nice cd-hit -d 0 -T 0 -M 0 -n $wordsize -c $identity -i $dir/proteins -o $dir/out";
msg("Running: $cmd");
system("$cmd 1> $dir/cdhit.stdout")==0 or err("Failed to run: $cmd");

#system("ls -lsa $dir");
#system("cp $dir/out.clstr .");

my %c;
my $cn=-1;
my $n=0;
open CLUSTERS, '<', "$dir/out.clstr";
while (<CLUSTERS>) {
  chomp;
  if (m/^>Cluster\s+(\d+)/) {
    $cn=$1;
    $n++;
  }
  elsif (m/^\d+\s+\d+aa, >(\d)~(.*?)\.\.\./) {
#    print "cn=$cn n=$n s=$1 d=$2\n";
    push @{ $c{$cn}{$1} }, $2;
  }
}
msg("Found $n clusters");
print Dumper(\%c) if $verbose;

my @uniq;

for my $i (keys %c) {
  if ( !defined $c{$i}{1} ) {    # unique to set 2
    push @{$uniq[2]}, @{ $c{$i}{2} }; # expand array (paralogs)
  }
  elsif ( !defined $c{$i}{2} ) { # unique to set 1
    push @{$uniq[1]}, @{ $c{$i}{1} }; # expand array (paralogs)
  }
}

for my $s (1,2) {
  printf "### Unique to %s - %d proteins\n", $set[$s], scalar(@{$uniq[$s]});
  my $count=0;
  for my $id ( @{$uniq[$s]} ) {
    printf "%d,%s,%s\n", ++$count, $id, $desc{$s}{$id};
  }
}


#>Cluster 0
#0       10421aa, >1~cds01403... at 95.06%
#1       10498aa, >2~gi|57285565|gb|AAW37659.1|... *
#>Cluster 1
#0       937aa, >1~cds02196... at 99.68%
#1       1510aa, >1~cds02197... at 99.74%
#2       2478aa, >2~gi|57286364|gb|AAW38458.1|... *


#----------------------------------------------------------------------

sub msg { print STDERR @_,"\n"; }
sub err { msg("ERROR: ", @_); exit -1; }

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"identity=f",  VAR=>\$identity, DEFAULT=>0.75, DESC=>"Percent identity for clustering"},
    {OPT=>"1=s",  VAR=>\$set[1], DEFAULT=>'', DESC=>"Protein set 1 [FASTA]"},
    {OPT=>"2=s",  VAR=>\$set[2], DEFAULT=>'', DESC=>"Protein set 2 [FASTA]"},
  );

  @ARGV or usage();

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] -1 alpha.faa -2 beta.faa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
