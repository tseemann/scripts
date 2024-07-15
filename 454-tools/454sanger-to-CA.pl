#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
#use Bio::Seq::Quality;

#>pdd4apcr001a04.b1 template=pdd4apcr001a04 dir=fwd trim=46-692
#gggaaattctgtaatgcggcggatttcgtaactgcctgttgaggtaagaaaaggccagct
#>tphlaa01p1001.sl1.b1 template=tphlaa01p1001 dir=fwd trim=13-764 library=tphl
#tgaatgatatacagcttgaattcgttacgaattctctagatatcgctcaatactgaccat
#>tphlaa01p1001.sr2.g1 template=tphlaa01p1001 dir=rev trim=1-730 library=tphl
#gagcgatatctagagaattcgtaacgaattcaagcttgatatcattcaggacgagcctca

my(@Options, $verbose);
setOptions();

for my $file (@ARGV) {
  print STDERR "Reading: $file\n";
  print STDERR "Writing: $file.vectorclear\n";
  open my $VEC, '>', "$file.vectorclear";
  my $mates=0;
  my $reads=0;
  my %mate;
  my $in = Bio::SeqIO->new(-file=>$file, -format=>'fasta');
  while (my $seq = $in->next_seq) {
    $reads++;
    if ($seq->desc =~ m/\btrim=(\d+)-(\d+)/) {
      print $VEC $seq->display_id." $1 $2\n";
    }
    if ($seq->desc =~ m/\btemplate=(\S+)/) {
      push @{$mate{$1}}, $seq->display_id;
    }
  }
  close $VEC;
  print STDERR Dumper(\%mate) if $verbose;

  print STDERR "Writing: $file.matepairs\n";
  open my $MATES, '>', "$file.matepairs";
  for my $id (keys %mate) {
    my @pair = @{$mate{$id}};
    next unless @pair == 2;
    print $MATES "@pair\n";
    $mates++;
  }
  close $MATES;
  printf STDERR "$file - found $mates mates from $reads reads, %d unmated.\n",
    ($reads - 2*$mates);
  
}



#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
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
  print "Usage: $0 [options] 454Sanger.fna\n";
  print "Output: 454Sanger.fna.{matepairs,vectorclear} for convert-fasta-to-v2.pl\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
