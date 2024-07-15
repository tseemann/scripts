#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use List::Util qw(sum);

my(@Options, $verbose, $ftype);
setOptions();

#tRNA            655..730
#                /gene="tRNA-Leu(UUR)"
#                /anticodon=(pos:678..680,aa:Leu)
#                /product="transfer RNA-Leu(UUR)"

my $gbk = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'genbank');

while (my $seq = $gbk->next_seq) 
{
  my %trna;
  my $count=0;
  for my $f (grep { $_->primary_tag eq 'tRNA' } $seq->get_SeqFeatures)  {
    my $p = TAG($f, 'product').' '.TAG($f, 'gene');
    print "$p\n";
    my $anti = '';
    if ($p =~ m/RNA-\w{3,4}([AGTCU]{3})/i) {
      $anti = $1;
    }
    elsif ($p =~ m/\b([AGTCU]{3})\b/) {
      $anti = $1;
    }
    elsif (TAG($f,'anticodon') =~ m/pos:(\d+)\.\.(\d+)/) {
      $anti = $f->seq->subseq($1,$2);
    }    
    $anti ||= 'NNN';
    $anti = uc($anti);
    $anti =~ s/U/T/g; # convert to DNA
    $anti =~ tr/ATGC/TACG/; # convert anticodon to codon
    $trna{$anti}++;    
  }
  print STDERR Dumper(\%trna);
  
  my $bad = $trna{NNN} || 0;
  printf STDERR "Parsed '%s' sequence.\n", $seq->display_id;
  printf STDERR "Found %d tRNA genes.\n", sum(values %trna);
  printf STDERR "WARNING: Could not determine anti-codon for %d tRNA genes.\n", $bad if $bad;
  printf STDERR "Found %d distinct codons.\n", (scalar keys %trna)-$bad;
  printf STDERR "Lacking %d codons (excluding stop)\n", 63-(scalar keys %trna)+$bad;
  printf STDERR "Duplicated codons: %s\n", 
    join(' ', map { "$_($trna{$_})" } grep { $trna{$_} > 1 } keys %trna);

  my %use;

  for my $f (grep { $_->primary_tag eq $ftype } $seq->get_SeqFeatures) {
    my $lt = TAG($f, 'locus_tag');
    my $dna = $f->seq->seq;
    my @codon = unpack '(A3)*', $dna;
#    print STDERR "$dna\n@codon\n"; exit;
    map { $use{$_}++ } @codon;
    my @bad = grep { not $trna{$_} } @codon;
    print STDERR "$lt\n$dna\n@codon\n@bad\n" if @bad; 
    exit;
    print STDERR "\r$lt";
  }
  print "\r";
  print STDERR Dumper(\%use);
}

#----------------------------------------------------------------------

sub TAG {
  my($seqfeat, $tag) = @_;
  return '' unless $seqfeat->has_tag($tag);
  return join(' ~~~ ', $seqfeat->get_tag_values($tag));
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"f|ftype=s",  VAR=>\$ftype, DEFAULT=>'CDS', DESC=>"Verify codons in these feature types"},
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
  print "Synopsis: Checks that there are tRNAs for all codons in genes\n";
  print "Usage: $0 [options] sequence.gbk\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
