#!/usr/bin/env perl
use strict;
use Data::Dumper;
use IO::File;
use Fatal;
use List::Util qw(min max);
use FindBin;
use lib "$FindBin::Bin/";
use FastQ qw(read_fastq write_fastq write_fasta assert);

my(@Options, $verbose, $output, $muxstart, $muxend, $readstart, $proportion);
setOptions();

my $fname = shift @ARGV;
printf STDERR "Determining barcodes: $fname\n";

my %tag;
my $nread=0;
my $nwrote=0;
open my $rfh, '<', $fname;

while (not eof $rfh) {			# perldoc -f eof
  my $r = read_fastq($rfh);
  $nread++;
  print STDERR "\rRead: $nread" if $nread % 17951 == 0;
#  my $L = length $r->[1];
#  next unless $L > $readstart;
  my $idx = substr($r->[1], $muxstart-1, $muxend-$muxstart+1);
  $tag{$idx}++;
}
print STDERR "Read $nread tags.\n";

my $max = max(values %tag);
printf STDERR "Most popular tag occurred $max times. Accepting those > %d\n", $proportion*$max;
for my $t (keys %tag) {
  delete $tag{$t} unless $tag{$t} > $proportion*$max;
}
print Dumper(\%tag);

for my $t (keys %tag) {
  my $outf = sprintf $output, $t;
  print STDERR "Creating $outf for tag $t\n";
  open my $fh, '>', $outf;
  $tag{$t} = $fh;
}

my $rej = sprintf( $output, 'N'x($muxend-$muxstart+1) );
print STDERR "Creating $rej file for rejected reads\n";
open my $rejfh, '>', $rej;

print STDERR "Rewinding filehandle: $fname\n";
seek $rfh, 0, 0;
$nread = $nwrote = 0;

#my $writer = $fasta ? \*write_fasta : \*write_fastq;

while (not eof $rfh) {			# perldoc -f eof
  my $r = read_fastq($rfh);
  $nread++;
  print STDERR "\rRead: $nread" if $nread % 17951 == 0;
  my $L = length $r->[1];
  next unless $L > $readstart;
  my $idx = substr($r->[1], $muxstart-1, $muxend-$muxstart+1);
  $r->[1] = substr($r->[1], $readstart-1);
  $r->[2] = substr($r->[2], $readstart-1);
  write_fastq($tag{$idx}||$rejfh, $r);
  $nwrote++;
}
print STDERR "\n";

printf STDERR "Parsed %d reads. Discarded %d. Wrote %d (%.2f%%)\n", 
  $nread, 
  $nread-$nwrote, 
  $nwrote, ($nwrote*100/$nread);

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"o|output=s", VAR=>\$output, DEFAULT=>'s_%s.fq', DESC=>"Output template: %s=barcode"},
#    {OPT=>"f|fasta!", VAR=>\$fasta, DEFAULT=>0, DESC=>"FASTA, not FASTQ output"},
    {OPT=>"b|muxstart=i", VAR=>\$muxstart, DEFAULT=>1, DESC=>"Start base of mux tag"},
    {OPT=>"c|muxend=i", VAR=>\$muxend, DEFAULT=>4, DESC=>"End base of mux tag"},
    {OPT=>"s|readstart=i", VAR=>\$readstart, DEFAULT=>8, DESC=>"Start base of read"},
    {OPT=>"p|proportion=f", VAR=>\$proportion, DEFAULT=>0.25, 
      DESC=>"Minimum tag frequency to accept as proportion of maximal occuring tag"},
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
  print "Usage: $0 --output 'demux_%b.fq' pooled.fq \n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
