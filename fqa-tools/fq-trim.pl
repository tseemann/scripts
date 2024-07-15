#!/usr/bin/env perl
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fastq trim_fastq assert write_fasta);

my(@Options, $verbose, $fasta, $minq, $minlen, $maxpoly);
setOptions();

assert($minq >= -5 and $minq <= 40, "--minq must be between -5 and 40");

my $fname = "@ARGV" || '(stdin)';
print STDERR "Trimming $fname to Q$minq. Please be patient...\n";

my($nread,$ntrim,$nambig,$nshort,$npass,$npoly,$oldbp,$newbp) = (0)x8;

my $writer_func = $fasta ? \&write_fasta : \&write_fastq;

while (not eof () ) {         # please read "perldoc -f eof" !
  my $s = read_fastq(\*ARGV);
  $nread++;
  $oldbp += length($s->[1]);
  
  $s = trim_fastq($s, $minq);
  next unless defined $s;
  $ntrim++;

  my $L = length($s->[1]);
  
  next if $L < $minlen;
  $nshort++;
  
  next if $L > $maxpoly and $s->[1] =~ m/^(.)\1+$/;
  $npoly++;
    
  next if index($s->[1],'N') > 0;
  $nambig++;  
  
  $writer_func->(\*STDOUT, $s);
  $npass++;
  $newbp += $L; # length($s->[1]);
}
if ($nread) {
  printf STDERR "%d reads, trimmed(Q$minq) %d, unambig %d, size(${minlen}bp) %d, nopoly %d, passed %d (%.2f%%)\n",
    $nread, $ntrim, $nambig, $nshort, $npass, $npoly, $npass*100/$nread;
  printf STDERR "read %d bp, wrote %d bp (%.2f%%)\n",
    $oldbp, $newbp, $newbp*100/$oldbp;
}
else {
  print STDERR "No reads found in $fname.\n";
}    


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
    {OPT=>"minq=i", VAR=>\$minq, DEFAULT=>20, DESC=>"Minimum quality: -5 to 40"},
    {OPT=>"minlen=i", VAR=>\$minlen, DEFAULT=>16, DESC=>"Minimum read length in output"},
    {OPT=>"maxpoly=i", VAR=>\$maxpoly, DEFAULT=>8, DESC=>"Maximum homopolymer read"},
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
  print "Usage: $0 [options] reads.fq\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
