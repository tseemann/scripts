#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $input, $list, $true, $false);
setOptions();

my %is_true;

open my $listf, "<$list" or die "bad --list";
while (<$listf>) {
  chomp;
  next unless $_;
  $is_true{$_}++;
}
my $nlist = scalar keys %is_true;
print STDERR "$list contains $nlist IDs.\n";


my $inf    = Bio::SeqIO->new(-file=>"<$input", -format=>'Fasta') or die "bad --input";
my $truef  = Bio::SeqIO->new(-file=>">$true",  -format=>'Fasta') or die "bad --true";
my $falsef = Bio::SeqIO->new(-file=>">$false", -format=>'Fasta') or die "bad --false";

my($nin,$ntrue,$nfalse) = (0,0,0);

while (my $seq = $inf->next_seq) {
  $nin++;
  if ( $is_true{ $seq->display_id } ) {
    print STDERR "Found match: ",$seq->display_id, " ", $seq->description, "\n" if $verbose;
    $truef->write_seq($seq);
    $ntrue++;
  }
  else {
    $falsef->write_seq($seq);
    $nfalse++;
  }
} 

print STDERR "Read $nin, positives $ntrue (of $nlist), negatives $nfalse\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"h|help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"i|input=s",  VAR=>\$input, DEFAULT=>'', DESC=>"Input FASTA file"},
    {OPT=>"l|list=s",  VAR=>\$list, DEFAULT=>'', DESC=>"Input TXT list of IDs"},
    {OPT=>"t|true=s",  VAR=>\$true, DEFAULT=>'/dev/null', DESC=>"Output for matching IDs"},
    {OPT=>"f|false=s",  VAR=>\$false, DEFAULT=>'/dev/null', DESC=>"Output for mis-matching IDs"},
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
  print "Usage: $0 [options]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
