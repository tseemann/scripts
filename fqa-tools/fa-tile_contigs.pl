#!/usr/bin/env perl
use warnings;
use strict;
use File::Temp qw(tempdir);
use Bio::SeqIO;
use File::Spec;
use List::Util qw(sum);
use Fatal;

print STDERR "Running: $0 @ARGV\n";

my(@Options, $verbose, $distant, $nodesc, $ref, $contigs, $out);
setOptions();

my $aligner = $distant ? 'promer' : 'nucmer';
$ref && -r $ref or die "Please specify a readable --reference file";
$contigs && -r $contigs or die "Please specify a readable --contigs file";

for my $tool (qw(nucmer promer mummerplot gnuplot show-snps show-coords delta-filter)) {
  my($path) = qx(which $tool);
  $path or die "Can not find '$tool'. Please install the MUMmer package.";
  print STDERR "Found $tool: $path";
}

my $dir = tempdir(CLEANUP=>1);
#my $dir = 'tmp'; system("mkdir $dir");
print STDERR "Using temporary folder: $dir\n";

my $refseq = load_fasta_to_hash($ref);
scalar keys %$refseq > 0 or die "--reference contains no sequences!";

my $ctgseq = load_fasta_to_hash($contigs);
scalar keys %$ctgseq > 0 or die "--contigs contains no sequences!";

my $nctg = scalar keys %$ctgseq;
my $ctgbp = sum( map { $_->length } values %$ctgseq );

print STDERR "Aligning contigs to reference with $aligner (please be patient)\n";
run_cmd("$aligner -maxmatch -p '$dir/out' \Q$ref\E \Q$contigs\E".
        " 2> '$dir/$aligner.err'");
print STDERR "Finding optimal tiling...\n";
# the -u option is essential for uniqifying alignments to plasmids etc
run_cmd("delta-filter -r -q -u 20 '$dir/out.delta' > '$dir/out.filter'");
run_cmd("mummerplot -p '$dir/plot' --layout --png '$dir/out.filter'".
        " 2> '$dir/mummerplot.err' 1> '$dir/mummerplot.out'");
run_cmd("show-coords -r -B '$dir/plot.filter' > '$dir/plot.filter.coords'");

# record which contigs were appended to end rather than truly aligned
my %aligned;
open COORDS, "$dir/plot.filter.coords";
while (<COORDS>) {
  my @x = split m/\t/;
  $aligned{ $x[0] } = $x[5];  # what it aligned to
}
my $naligned = scalar keys %aligned;
my $alignedbp = sum( map { $ctgseq->{$_}->length } keys %aligned );
printf STDERR "Placed %d of %d bp from contigs (%.2f%%)\n", $alignedbp, $ctgbp, 100*$alignedbp/$ctgbp;
printf STDERR "Placed %d of %d contigs (%.2f%%)\n", $naligned, $nctg, 100*$naligned/$nctg;

# now use the mummerplot tiling (but which appended unplaced contigs)
my $fout = Bio::SeqIO->new(-file=>">$out", -format=>'fasta');
my $parsenow=0;
open GNUPLOT, "$dir/plot.gp";
while (<GNUPLOT>) {
  if (m/set ytics/) {
    $parsenow=1;
    next;
  }
  elsif ($parsenow and m/^\s+\"(\*)?([^\"]+)\"\s+(\d+),\s+\\/) {
    my($strand,$id,$pos) = ($1 ? '-':'+' , $2, $3);
    if (exists $ctgseq->{$id}) {    
      $ctgseq->{$id} = $ctgseq->{$id}->revcom if $strand eq '-';
      if (not $nodesc) {        
        $ctgseq->{$id}->desc($aligned{$id} ? "aligned $strand $aligned{$id}" : "unplaced");
      }
      $fout->write_seq( $ctgseq->{$id} );
      delete $ctgseq->{$id};
    }
  }
}

printf STDERR "Writing %d orphan contigs.\n", scalar keys %$ctgseq;
for my $orphan (sort keys %$ctgseq) {
  $ctgseq->{$orphan}->desc("orphan") unless $nodesc;
  $fout->write_seq( $ctgseq->{$orphan} );
}

print STDERR "Results in file: $out\n";
print STDERR "Done.\n";

#----------------------------------------------------------------------
# Useful functions

sub run_cmd {
  my($cmd) = @_;
  print STDERR "Running: $cmd\n";
  system($cmd)==0 or die "ERROR: could not run command '$cmd'";
}

sub load_fasta_to_hash {
  my($fname) = @_;
  my %hash;
  my $counter=0;
  my $bp=0;
  print STDERR "Reading FASTA: $fname\n";
  my $in = Bio::SeqIO->new(-file=>$fname, -format=>'fasta');
  while (my $seq = $in->next_seq) {
    my $id = $seq->id;
    $counter++;
    if (exists $hash{$id}) {
      print STDERR "WARNING: $fname has duplicate sequence ID: $id - appending '.$counter'\n";
      $id .= ".$counter";
    }
    $hash{$id} = $seq;
    $bp += $seq->length;
  }
  printf STDERR "Loaded %d bp in %d contigs\n", $bp, $counter;
  return \%hash;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,            DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"distant!",  VAR=>\$distant, DEFAULT=>0, DESC=>"Genomes are divergent at DNA level"},
    {OPT=>"nodesc!",  VAR=>\$nodesc, DEFAULT=>0, DESC=>"Don't put alignment details to FASTA description"},
    {OPT=>"ref=s",  VAR=>\$ref, DEFAULT=>'', DESC=>"Reference genome [FASTA]"},
    {OPT=>"contigs=s",  VAR=>\$contigs, DEFAULT=>'', DESC=>"Draft contigs [FASTA]"},
    {OPT=>"outfile=s",  VAR=>\$out, DEFAULT=>'/dev/stdout', DESC=>"Output file [FASTA]"},
  );

#  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  my(undef,undef,$exe) = File::Spec->splitpath($0);
  print "Synopsis:\n  Order and orient draft contigs against a reference genome\n";
  print "Usage:\n  $exe [options] -r ref.fa -c contigs.fa -o tiled.fa\n";
  print "Options:\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
