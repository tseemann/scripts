#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use List::MoreUtils qw(pairwise each_array);
use List::Util qw(sum);

no warnings qw(once); #Name "main::a" used only once: possible typo

my(@Options, $verbose, $prefix, $ambig, $mincov);
setOptions();

my $refname = "$ARGV[0]/reference.fa";
print STDERR "Loading: $refname";
my $refio = Bio::SeqIO->new(-file=>$refname, -format=>'fasta');
my $ref = $refio->next_seq->seq;
my $len = length($ref);
print STDERR " - $len bp (reference genome)\n";

# initialise with 1s - we will multiply by coverage at each point
# any zero coverage will zero it forever
my @cov = map { 1 } (1 .. $len);
my @totcov = map { 0 } (1.. $len);
my $numdir = 0;
my %snp;

print STDERR "Writing: $prefix.txt\n";
open my $logfh, '>', "$prefix.txt";
print $logfh tsv(qw(Genome MeanCov ZeroCov));

for my $ndir (@ARGV) {
  next unless -d $ndir;
  my($fn) = glob("$ndir/*-ambiguous-depth.userplot");
  next unless $fn;
  print STDERR "Loading: $fn - ";
  if (not -r $fn) {
    print STDERR " could not read, skipping.\n";
  }
  else {
    open my $fh, '<', $fn;
    my @depth = <$fh>;
    close $fh;
    printf STDERR "%d lines\n", scalar @depth;
    if (@depth != $len) {
      print STDERR "ERROR: $fn does not have length $len bp - skipping\n";
      next;
    }
    # the trick! zero out any position that has no coverage
    # where coverage is defined as "less than $mincov".
    chomp @depth;
    @cov = pairwise { $a * ($b > $mincov) } @cov, @depth;
    # keep running total for colouring non-core
    @totcov = pairwise { $a + ($b > $mincov) } @totcov, @depth;
    print STDERR "Processed: $ndir - core is now ",sum(@cov),"/$len\n";
    my $mcov = sprintf "%.2f", sum(@depth) / scalar(@depth);
    my $zcov = sprintf "%.2f", scalar(grep { $_ < 1 } @depth) / scalar(@depth);
    printf $logfh tsv($ndir, $mcov, $zcov);
    load_snps($ndir);
    $numdir++;
  }
}

print STDERR "Writing: $prefix.snps.csv\n";
open my $snpfh, '>', "$prefix.snps.csv";

print STDERR "Writing: $prefix.masked.fasta\n";
open my $seqfh, '>', "$prefix.masked.fasta";

print STDERR "Writing: $prefix.coreness.userplot\n";
open my $plotfh, '>', "$prefix.coreness.userplot";

printf $seqfh ">$prefix core genome $len bp\n";
my $ncore=0;
my @gid = (sort keys %snp);
print $snpfh tsv('Position', 'Reference', @gid);
my $snprows=0;

my %tree;
for my $i (0 .. $len-1) {
  print STDERR "\rProcessing: $i/$len bp" if not $verbose and $i%3971==0;
  # the userplot will have #genomes mismatched to make non-core
  print $plotfh ($numdir-$totcov[$i]),"\n";
  # get base call
  my $base = substr($ref,$i,1);
  if ($cov[$i]) {
    print $seqfh uc($base);
    $ncore++;
    my $SAME = $base;
#    my $SAME = ' ';
    my @snp = map { $snp{$_}{$i+1} || $SAME } @gid;
    if (grep { $_ ne $SAME } @snp) {
      print $snpfh tsv($i+1, $base, @snp);
      $tree{'Reference'} .= $base;
      my $ea = each_array(@gid, @snp);
      while ( my ($g, $s) = $ea->() ) {
#        $tree{$g} .= $s;
        $tree{$g} .= substr $s, 0, 1; # hack to handle indels > 1bp
      }
      $snprows++;
    }
  }
  else {
    print $seqfh ($ambig ? 'N' : lc($base));
  }
  print $seqfh "\n" if (($i+1)%60==0);
  
}
print $seqfh "\n";

printf STDERR "\nCore genome is $ncore/$len bases (%.1f%%)\n", $ncore*100/$len;
printf $logfh "\nCore genome is $ncore/$len bases (%.1f%%)\n", $ncore*100/$len;

print STDERR "Writing: $prefix.tree.fasta\n";
#open my $treefh, '>', "$prefix.tree.fasta";
my $treeio = Bio::SeqIO->new(-file=>">$prefix.tree.fasta", -format=>'fasta');
my $desc = sprintf "core genome %d SNPs across %d genomes", $snprows, 1+scalar(@gid);
for my $g (sort keys %tree) {  # should be same as (@gid) 
#  print $treefh ">$g\n",$tree{$g},"\n";
  $treeio->write_seq(Bio::Seq->new(
    -id=>$g, 
    -seq=>$tree{$g}, 
    -desc=>$desc,
  ));
}

#print STDERR Dumper(\%snp);

#----------------------------------------------------------------------

sub tsv {
  join("\t", @_)."\n";
}

#----------------------------------------------------------------------
#Sequence        Position in reference   Change type     Old     New     Evidence
#Saa_6008        10866   substitution    T       C       "C"x39

sub load_snps {
  my($dir) = @_;
  my $repfn = "$dir/report.txt";
  -r $repfn or die "ERROR: could not read '$repfn'";
  print STDERR "Loading: $repfn";
  open my $fh, '<', $repfn;
  while (<$fh>) {
    my @x = split m/\t/;
    next unless @x==6 and $x[1] =~ m/^\d+$/;
    my $change = ($x[2] =~ m/insert/) ? lc($x[4]) : uc($x[4]);
    $snp{$dir}{$x[1]} = $change;
  }
  print STDERR " - read ", scalar(keys %{$snp{$dir}}), " SNPs\n";
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"p|prefix=s",  VAR=>\$prefix, DEFAULT=>'uniqueome', DESC=>"Prefix for output files: .fna .csv .txt"},
    {OPT=>"n|ambig!",  VAR=>\$ambig, DEFAULT=>0, DESC=>"Ouput Ns rather than lowercase AGTC for non-core"},
    {OPT=>"m|mincov=i",  VAR=>\$mincov, DEFAULT=>0, DESC=>"Threshold to count as 'no coverage'"},
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
  print "Usage: $0 [options] nesoniDir1 [nesoniDir2 ...]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
