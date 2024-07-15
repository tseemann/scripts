#!/usr/bin/perl -w
use strict;
#use threads;
use Data::Dumper;
use List::Util qw(min max);
use List::MoreUtils qw(uniq);

my(@Options, $verbose, $force, $cpu, $makefile,
             $k_patt, $mincov_patt, $expcov_patt, $maxcov_patt, $insert_patt,
	     $opts);
setOptions();

my $dir = shift @ARGV or die "need velveth_dir";
-r "$dir/Sequences" or die "can't see $dir/Sequences";
print STDERR "Using velveth result: $dir\n";

#$cpu = num_cpu() if $cpu <= 0;
#print STDERR "Found $cpu CPUs. (only using 1)\n";

if ($insert_patt) {
  my($mu,$sd) = decode_pair($insert_patt);
  $opts .= " -ins_length $mu" if $mu;
  $opts .= " -ins_length_sd $sd" if $sd;
#  $opts .= " -scaffolding yes";
# THESE WILL CAUSE BUGS!
#  $opts .= " -ins_length $mu -ins_length2 $mu" if $mu;
#  $opts .= " -ins_length_sd $sd -ins_length_sd2 $sd" if $sd;
}
print STDERR "Global velveth parameters: $opts\n";

runcmd("velvetg $dir $opts > /dev/null") unless -r "$dir/stats.txt";
#runcmd("sort -k 2 -n -r $dir/stats.txt | head -10 | cut -f2,7");

my $est_cov = estimate_expcov("$dir/stats.txt");
print STDERR "Estimated k-mer coverage -exp_cov: $est_cov\n";

#my $est_bases = estimate_bases("$dir/Sequences");
#print STDERR "Estimated raw depth: $est_depth\n";

my @expcov = $expcov_patt ? decode_range($expcov_patt) : ('auto');
#                          : uniq($est_cov, int(0.7*$est_cov), int(0.85*$est_cov), int(1.15*$est_cov), int(1.3*$est_cov));
print STDERR "expcov: @expcov\n";

my @mincov = $mincov_patt ? decode_range($mincov_patt) : ('auto');
#                          : uniq(0 .. int(0.3 * $est_cov));
print STDERR "mincov: @mincov\n";

if ($makefile) {
  generate_makefile();
}

#my @maxcov = $maxcov_patt ? decode_range($maxcov_patt) 
#                          : uniq(100 * $est_cov);
#print STDERR "maxcov: @maxcov\n";

for my $expc (@expcov) {
  for my $minc (@mincov) {
    my $suffix = "exp${expc}_min${minc}";
    print STDERR "Processing: $dir | $suffix\n";
    runcmd("velvetg $dir -cov_cutoff $minc -exp_cov $expc $opts 1> $dir/stdout.$suffix");    
    runcmd("cp -f $dir/contigs.fa $dir/$suffix.fna");
    runcmd("cp -f $dir/stats.txt $dir/stats.$suffix");
    runcmd("fa $dir/$suffix.fna");
  }
}

runcmd("fa -e $dir/*.fna");


#  for my $id (1 .. $fileno) {
#    my @cmd = ('nice','phmmer', '-E', $evalue, '--notextw', '-o', "$hitfile.$id", "$id.fasta", $seqfile);
#    print STDERR "Spawning thread $id: @cmd\n";
#    my $thr = threads->new(sub { system(@_) }, @cmd);
#  }
#  for my $thr (threads->list) {
#    print STDERR "Waiting for thread ",$thr->tid," to complete.\n";
#    $thr->join;
#  }

sub runcmd {
  my(@cmd) = @_;
  print STDERR "Running: @cmd\n";
  system(@cmd)==0 or die "ERROR: $?";
}

sub decode_range {
  my($patt) = @_;
  my @p = split m/\s*,\s*/, $patt;
  my @list;
  for my $p (@p) {
    $p =~ m/ ^ \s* (\d+) (?:-(\d+))? (?:\/(\d+))? \s* $ /xms;
    my $start = defined($1) ? $1 : die "need FROM in $p";
    my $end = $2 || $start;
    my $step = $3 || 1;
    if ($end < $start) { ($end,$start) = ($start,$end) }
    for (my $i=$start; $i <= $end; $i+=$step) {
      push @list, $i;
    }
  }
  return @list;
}

sub decode_pair {
  my($p) = @_;
  $p =~ m/(\d+)(?::(\d+))?/;
  return ($1,$2);
}

sub num_cpu {
  my($n) = qx(grep -c ^bogomip /proc/cpuinfo);
  chomp $n;
  return $n;
}

sub estimate_expcov {
  my($statfile) = @_;
  my($mode1) = qx(cat $statfile | sort -k 2 -n -r | cut -f 7 | math numeric | math quantize 1 | math above 6 | math mode);
  chomp $mode1;
  my($mode2) = qx(cat $statfile | sort -k 2 -n -r | cut -f 9 | math numeric | math quantize 1 | math above 6 | math mode);
  chomp $mode2;
  return $mode1 + $mode2;  # add -short1 and -short2 coverages
}  

sub estimate_bases {
  my($seqfile) = @_;
  my($wc) = qx(grep -v '^>' $seqfile | wc);
  chomp $wc;  
  my($lines,undef,$chars) = split ' ', $wc;
  return $chars-$lines;
}  

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"opts=s",  VAR=>\$opts, DEFAULT=>'', DESC=>"common velvetg -options. Type 'velvetg' to list them."},  
#    {OPT=>"f|force!",  VAR=>\$force, DEFAULT=>0, DESC=>"Force over-write of existing outputs"},
#    {OPT=>"cpu=i",  VAR=>\$cpu, DEFAULT=>8, DESC=>"Number of CPUs to use (NOT IMPLEMENTED)"},  
#    {OPT=>"k=s",  VAR=>\$k_patt, DEFAULT=>'29-33/2', DESC=>"velveth kmer value: FROM-TO/STEP,FROM-TO,VALUE"},  
    {OPT=>"cov_cutoff=s",  VAR=>\$mincov_patt, DEFAULT=>0, DESC=>"-cov_cutoff FROM-TO/STEP,FROM-TO,VALUE (0=auto)"},  
    {OPT=>"exp_cov=s",  VAR=>\$expcov_patt, DEFAULT=>0, DESC=>"-exp_cov FROM-TO/STEP,FROM-TO,VALUE (0=auto)"},  
#    {OPT=>"maxcov=s",  VAR=>\$maxcov_patt, DEFAULT=>0, DESC=>"-max_coverage FROM-TO/STEP (0=auto)"},  
    {OPT=>"insert=s",  VAR=>\$insert_patt, DEFAULT=>'250:45', DESC=>"-ins_length{_sd} MEAN:SD (0=auto)"},  
#    {OPT=>"insert=s",  VAR=>\$insert_patt, DEFAULT=>'250:45', DESC=>"-ins_length{_sd} MEAN:SD (0=auto)"},  
#    {OPT=>"insert=s",  VAR=>\$insert_patt, DEFAULT=>'250:45', DESC=>"-ins_length{_sd} MEAN:SD (0=auto)"},  
#    {OPT=>"makefile!",  VAR=>\$makefile, DEFAULT=>0, DESC=>"Output Makefile to stdout"},  
  );

#  @ARGV > 1 or usage();

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();


  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] velveth_dir/ \n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
