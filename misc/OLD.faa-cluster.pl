#!/usr/bin/env perl
use strict;
use threads;
use Bio::SeqIO;
#use File::Spec;
use Regexp::Common;
use Data::Dumper;
use List::Util qw(min max);
use List::PowerSet qw(powerset);
use List::MoreUtils qw(uniq);
#use sort 'stable';
use Text::CSV;

my(@Options, $verbose, $force, $cpu, $inflation, $mclexe, $evalue, $prefix);
setOptions();
$|=1;

my $SEP = '.'; # DO NOT CHANGE THIS OR int() WILL FAIL

my $dir = shift @ARGV;
if (-e $dir) {
  print STDERR "$dir already exists.\n";
}
else {
  mkdir $dir or die "could not mkdir '$dir'";
}

my($numseq) = qx(cat @ARGV | grep -c '^>');
chomp $numseq;
print STDERR "Estimated number of sequences: $numseq\n";

my $threshold = max(1, int($numseq / $cpu));
print STDERR "Splitting sequences into $cpu batches of ~$threshold sequences.\n";
my $fileno=0;
my $filecount=0;

my $seqfile = "all.fasta";
open my $FOUT, '>', "$dir/$seqfile";

my $splitfh;
my @tempfafiles;

my %seq;
my $genome=0;
for my $file (@ARGV) {
  my $protein=0;
  print STDERR "Loading: $file\n";
  my $in = Bio::SeqIO->new(-file=>$file, -format=>'Fasta');
  $file =~ m/^(\w+)/;
  my $prefix_str = $1 || $file;
  $genome++;
  while (my $seq = $in->next_seq) {
    my $s = uc $seq->seq; # uppercase
    $s =~ s/\*$//;        # remove final stop codon
    $s =~ s/\*/X/g;       # repleace in-frame stops with X
    $seq->display_id( $prefix_str . ':' . $seq->display_id ) if $prefix;
    $protein++;
    my $id = sprintf "%d%s%06d", $genome, $SEP, $protein;
    $seq{$id} = $seq;
    my $data = sprintf ">$id ".$seq->desc."\n$s\n", $genome, $protein;
    printf $FOUT $data;
    if (++$filecount > $threshold or $fileno==0) {
      $fileno++;
      $filecount=0;
      my $fname = "$fileno.fasta";
      print STDERR "Writing: $fname\n";
      open $splitfh, '>', "$dir/$fname";
      push @tempfafiles, $fname;
    }
    print $splitfh $data;
  } 
}
print STDERR "Stored ",scalar(keys %seq), " sequences.\n";

chdir $dir;
print STDERR "Changed directory: ", qx(pwd);

#my $hitfile = 'all.phmmer';
#my @cmd = ('nice', 'time', 'phmmer', '--notextw', '-o', $hitfile, $seqfile, $seqfile);
#print STDERR "Running: @cmd\n";
#system(@cmd);

my $hitfile = 'all.phmmer';
my $done_phmmer = !(!-r $hitfile or $force);

if (not $done_phmmer) {
  for my $id (1 .. $fileno) {
    my @cmd = ('nice','phmmer', '-E', $evalue, '--notextw', '-o', "$hitfile.$id", "$id.fasta", $seqfile);
    print STDERR "Spawning thread $id: @cmd\n";
    my $thr = threads->new(sub { system(@_) }, @cmd);
  }
  for my $thr (threads->list) {
    print STDERR "Waiting for thread ",$thr->tid," to complete.\n";
    $thr->join;
  }
  print STDERR "Combining hitfiles...\n";
  my @hitfiles = map { "$hitfile.$_" } (1 .. $fileno);
  system("cat @hitfiles > $hitfile");
  print STDERR "Deleting partial hit files: @hitfiles\n";
  unlink @hitfiles;
}
else {
  print STDERR "Re-using existing $hitfile ; use --force, or delete $hitfile, to re-run \n";
}  


print STDERR "Deleting partial fasta files: @tempfafiles\n";
unlink @tempfafiles;

print STDERR "Parsing: $hitfile\n";
open my $hitfh, '<', $hitfile; 
my $result = parse_hmmer3_file($hitfh);
#print Dumper($result);

sub canon {
  return join('~', sort { $a <=> $b } @_);
}

my %seen_pair;
my $mclinfile = 'mcl.in';
open my $mclin, '>', $mclinfile;
for my $q (@$result) {
  for my $h (@$q) {
    next if $seen_pair{ canon($h->{qid}, $h->{hid}) } ++;
    # print Dumper($h);
#    my $score = log( max(1, $h->{score}) );
    my $score = $h->{score};
    next unless $score > 0;
    print $mclin join(' ', $h->{qid}, $h->{hid}, $score), "\n";
  }
}

my $mcloutfile = 'mcl.out';
my @cmd = ($mclexe, $mclinfile, '-V', 'silent', '--abc', '-o', $mcloutfile, '-I', $inflation);
print STDERR "Running: @cmd\n";
system(@cmd);

print STDERR "Reading clusters ...\n";
my %is_alone = (map { ($_ => 1) } keys %seq); # single until proven not so
my @cluster;
open my $MCL, '<', $mcloutfile;
while (<$MCL>) {
  chomp;
  my @id = sort { $a <=> $b } split ' ';
  push @cluster, \@id;
  map { $is_alone{$_}=0 } @id;
#  print STDERR "$. @id\n";
}
print STDERR "Adding singletons ...\n";
for my $id (sort grep { $is_alone{$_} } keys %is_alone) {
#  print STDERR "singleton $id\n";
  push @cluster, [ $id ];
}
#print Dumper(\@cluster);

my @set = @{ powerset(1 .. $genome) };
@set = grep { scalar(@$_) > 0 } @set; # remove null group
@set = sort { scalar(@$a) <=> scalar(@$b) or join('',@$a) <=> join('',@$b) } @set;
#print STDERR Dumper(\@set);

sub cluster_in_set {
  my($set, $ids) = @_;
  my $setmark = join '', uniq @$set;
  my @ids = map { int } @$ids;
#  map { s/$SEP.*$//g } @ids;
  my $idmark = join '', uniq @ids;
#  print STDERR Dumper($set,$setmark,$ids,$idmark), "--------------\n";
#  print STDERR "$idmark == $setmark ?\n";
  return $idmark eq $setmark;
}

my $csv = Text::CSV->new;
my(@header,@body);
push @header, csv( qw(Venn NumClusters) );
my $bigcount=0;
for my $set (@set) {
  my $count=0;
  push @body, "\n", csv('SET', join('+',map { $ARGV[$_-1] } @$set) );
  push @body, csv('TOTAL','SET',map { $ARGV[$_-1] } (1 .. $genome));
  for my $cluster (@cluster) {
    next unless cluster_in_set($set, $cluster);
    $count++;
    $bigcount++;
    my @line = ($bigcount, $count, ("-")x$genome, "len", "desc");
    for my $n (@$set) {
      my @id = grep { m/^$n$SEP/ } @$cluster;
      $line[1+$n] = join(' ', map { $seq{$_}->display_id } @id);
      $line[1+$genome+1] = $seq{$id[0]}->length;
      $line[1+$genome+2] = $seq{$id[0]}->desc;
    }
    push @body, csv(@line);
  }
  push @body, "NONE" if $count==0;
  push @header, csv(join('+',@$set), $count);
}

open my $REPORT, '>', 'orthologs.csv';
print $REPORT @header, "\n";
print $REPORT @body;
close $REPORT;

print STDERR "Results in: $dir/\n";
exit;

sub csv {
  return unless @_;
  $csv->combine(@_) or die $csv->error_diag;
  $csv->string."\n";
}


#Query:       L550/2/891/139  [L=251]

#>> L550/2/891/139  ParA-like protein {LBL_4000} {SPN3821}
#  #    score bias c-Evalue i-Evalue hmmfrom hmm to        alifrom ali to     envfrom env to      acc
#---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#  1 !   42.0   0.0   2.2e-12   8.6e-10       2      83 ..     439     520 ..     437     521 .. 0.95
#  2 !   15.1   0.0   0.00057      0.22      14      82 ..     837     911 ..     827     914 .. 0.82
#  3 ?    3.6   0.0       2.4   9.2e+02       6      35 ..    1205    1235 ..    1202    1258 .. 0.81
#  4 !   22.9   0.0   2.1e-06   0.00081      13      79 ..    1312    1380 ..    1304    1385 .. 0.81
#  5 !   43.3   0.4   8.9e-13   3.4e-10       1      82 [.    1799    1888 ..    1799    1891 .. 0.89
#  6 !   19.2   0.0   3.1e-05     0.012       6      73 ..    1904    1966 ..    1900    1976 .. 0.90
#  7 !   12.0   0.0    0.0053       2.1       1      77 [.    1993    2098 ..    1993    2107 .. 0.74
#  0 1    2       3      4           5        6      7  8      9       10   11    12      13  14  15
# The ! or ? symbol indicates whether this domain does or does not satisfy both
# per-sequence and  per-domain inclusion thresholds.

# .. meaning both ends of the alignment ended internally, and 
# [] means both ends of the alignment were full-length flush to the ends of the query or target, and 
# [. and .] mean only the left or right end was flush/full length.


sub parse_hmmer3_file {
  my($fh) = @_;  
  local $/ = "//\n"; # input record separator
  my @res;
  while (<$fh>) {
    push @res, parse_hmmer3_result($_);
  }
  return \@res;
}

sub parse_hmmer3_result {
  my($s) = @_;
  my @res;
  $s =~ m/^Query:\s+(\S+)\s+\[L=(\d+)\]/xms or die "no Query: in: $s";
  my $qid = $1;
  my $qlen = $2;
  while ($s =~ m/ (^>> .*?)$ ^ $ /xmsg) {
#    print "BLOCK:\n$1\n";
    my $hit = parse_hmmer3_hit($1);
    next unless $hit;
    next if $hit->{hid} eq $qid;  # REMOVE SELF HITS
    $hit->{qid} = $qid;
    $hit->{qlen} = $qlen;
    push @res, $hit;    
  }
  return [ sort { $a->{evalue} <=> $b->{evalue} } @res ];
}

sub parse_hmmer3_hit {
  my($s) = @_;
  my @s = split m/\n/, $s;
  $s[0] =~ m/^>>\s+(\S+)/ or die "need >> ID in: $s[0]";
  my $id = $1;
#  die "no s(3) line in $s[0]" unless defined $s[3];
  return if not defined $s[3];  # no hits found
  my @h = split ' ', $s[3];
  die "need 16 cols in: $s[3]" unless @h == 16;
  my %hit = (
    hid    => $id,
    score  => $h[2],
    evalue => $h[5],
    qbegin => $h[6],
    qend   => $h[7],
    sbegin => $h[9],
    send   => $h[10],
  );
  return \%hit;
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"f|force!",  VAR=>\$force, DEFAULT=>0, DESC=>"Force over-write of existing outputs"},
    {OPT=>"I=f",  VAR=>\$inflation, DEFAULT=>1.5, DESC=>"MCL inflation value"},
    {OPT=>"e|evalue=f",  VAR=>\$evalue, DEFAULT=>1E-3, DESC=>"Similarity cut-off : evalue"},
    {OPT=>"mcl=s",  VAR=>\$mclexe, DEFAULT=>'mcl', DESC=>"mcl binary"},  
    {OPT=>"cpu=i",  VAR=>\$cpu, DEFAULT=>8, DESC=>"Number of CPUs to use"},  
    {OPT=>"prefix!",  VAR=>\$prefix, DEFAULT=>0, DESC=>"Prefix FASTA IDs with filename"},  
  );

  @ARGV >= 3 or usage();

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] OUTDIR/ gen1.faa gen2.faa [ gen3.faa ... ]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
