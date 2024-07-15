#!/usr/bin/env perl
use strict;
use SWISS::Entry;
use SWISS::KW;
#use SWISS::OC;
use Data::Dumper;

my(@Options, $verbose, $frag, $evlev, $minlen, $sep, $blank, $term);
setOptions();

# Read an entire record at a time
local $/ = "\n//\n";

my $in=0;
my $out=0;

while (<ARGV>) 
{
  # Read the entry
  my $entry = SWISS::Entry->fromText($_);
  $in++;

  # Immediately reject partial genes
  next if not $frag and $entry->isFragment;
  next if not $frag and $entry->DEs->hasFragment;

  # Too short?
  next if $minlen > 0 and length($entry->SQs->seq) < $minlen;
 
  # Reject on evidence level
  # grep ^PE uniprot_sprot.dat | sort | uniq -c
  #  74284 PE   1: Evidence at protein level;
  #  67762 PE   2: Evidence at transcript level;
  # 376894 PE   3: Inferred from homology;
  #  14424 PE   4: Predicted;
  #   1884 PE   5: Uncertain;
  if ($evlev < 5) { 
    $entry->PE->text =~ m/^(\d+)/;
    next unless $1 <= $evlev; 
  }

  # Only specified organism class
  if ($term) {
    my $tax = $entry->OCs->list or next;
    next unless grep { $_ eq $term } @$tax ;
  }

  # /gene code  
  my $gene = $entry->GNs->getFirst || '';
  $gene = '' if $gene =~ m/\d{2}/ or $gene =~ m/\./;

  my $ec = ''; 
  my $prod = 'hypothetical protein'; 

  if (1) {  
    for my $de ($entry->DEs->elements) {
      if ($de->type eq 'EC') {
	$ec = $de->text;
	$ec =~ s/^\D*//;
	last;
      }
      elsif ($de->type eq 'Full') {
	$prod = $de->text;
	if ($prod =~ m/^UPF\d|^Uncharacterized protein|^ORF|^Protein /) {
	  $prod = 'hypothetical protein';
	}
      }
    }
  }

  $ec ||= $blank;
  $gene ||= $blank;
  $prod ||= $blank;
  
  print STDERR join("\t", $entry->AC, $ec, $gene, $prod), "\n";
  print ">", $entry->AC, " $ec$sep$gene$sep$prod\n", $entry->SQs->seq, "\n";
  
  $out++;

}

#print STDERR "\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"separator=s",   VAR=>\$sep, DEFAULT=>'~~~', DESC=>"Separator for gene/EC/product"},
    {OPT=>"blank=s",   VAR=>\$blank, DEFAULT=>'', DESC=>"Replace empty gene/EC/product with this"},
    {OPT=>"evidence=i",   VAR=>\$evlev, DEFAULT=>6, DESC=>"1=prot 2=mrna 3=homol 4=pred 5=unsure"},
    {OPT=>"fragments!",   VAR=>\$frag, DEFAULT=>0, DESC=>"Include 'DE Flags: Fragment;' entries"},
    {OPT=>"minlen=i",   VAR=>\$minlen, DEFAULT=>0, DESC=>"Minimum peptide length"},
    {OPT=>"term=s",   VAR=>\$term, DEFAULT=>'Bacteria', DESC=>"Lineage must contain this term"},
  );

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] <uniprot.dat>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
