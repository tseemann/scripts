#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $fuzzy, $input, $list, $true, $false);
setOptions();

open my $listf, "<$list" or die "bad --list";
my @list = grep { m/./ } <$listf>;
my $nlist = scalar(@list);
print STDERR "$list contains $nlist patterns.\n" if $verbose;
chomp(@list);
close $listf;

my $pattern;
if ($fuzzy) {
  $pattern = join('|', map { quotemeta } @list);
}
else {
  $pattern = join('|', map { '^'.quotemeta($_).'$' } @list);
}
print STDERR "Pattern: $pattern\n" if $verbose;

my $inf    = Bio::SeqIO->new(-file=>"<$input", -format=>'Fasta') or die "bad --input";
my $truef  = Bio::SeqIO->new(-file=>">$true",  -format=>'Fasta') or die "bad --true";
my $falsef = Bio::SeqIO->new(-file=>">$false", -format=>'Fasta') or die "bad --false";

my($nin,$ntrue,$nfalse) = (0,0,0);

while (my $seq = $inf->next_seq) {
  $nin++;
  if ($seq->display_id =~ m/($pattern)/) {
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
    {OPT=>"z|fuzzy!",  VAR=>\$fuzzy, DEFAULT=>0, DESC=>"Allow fuzzy/partial matching IDs"},
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
