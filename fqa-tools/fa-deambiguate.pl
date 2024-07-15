#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $doc, $enns);
setOptions();

if ($doc) {
  print STDERR $_ for (<DATA>);
  exit;
}

my %is = (
  'm' => [ qw(a c) ],        #a or c
  'r' => [ qw(a g) ],        #a or g
  'w' => [ qw(a t) ],        #a or t
  's' => [ qw(c g) ],        #c or g
  'y' => [ qw(c t) ],        #c or t
  'k' => [ qw(g t) ],        #g or t
  'v' => [ qw(a c g) ],      #a or c or g; not t
  'h' => [ qw(a c t) ],      #a or c or t; not g
  'd' => [ qw(a g t) ],      #a or g or t; not c
  'b' => [ qw(c g t) ],      #c or g or t; not a
  'n' => [ qw(a c g t) ],    #a or c or g or t
);

sub deambig {
  my($c) = @_;
  return 'N' if $enns;
  $c = lc($c);
  die "bad letter '$c'" unless exists $is{$c};
  return $is{$c}->[ int(rand scalar( @{$is{$c}}) ) ];
} 

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $nchange=0;
my $nwrote=0;

while (my $seq = $in->next_seq) {
#  $out->write_seq($seq);
  if ($seq->alphabet eq 'dna') {
    my $s = $seq->seq;
    $s =~ s/([mrwsykvhdbn])/deambig($1)/ieg;
    $nchange++ if $seq->seq ne $s;
    $seq->seq($s);
  }
  $out->write_seq($seq);
  $nwrote++;
} 

print STDERR "Read $nwrote sequences, changed $nchange sequences.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"d|doc!",  VAR=>\$doc, DEFAULT=>0, DESC=>"Show documentation"},
    {OPT=>"n!",  VAR=>\$enns, DEFAULT=>0, DESC=>"Use Ns rather than random disambiguation"},
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
  print "Usage: $0 [options] ambig_dna.fasta > unambig_dna.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__
http://www.ncbi.nlm.nih.gov/collab/FT/

7.5.1 Nucleotide base codes (IUPAC)

Authority       Nomenclature Committee of the International Union of 
                Biochemistry 
Reference       Cornish-Bowden, A.  Nucl Acid Res 13, 3021-3030 (1985)
Contact         EMBL-EBI
Scope           Location descriptors 

Listing

        Symbol  Meaning
        ------  -------

        a       a; adenine
        c       c; cytosine
        g       g; guanine
        t       t; thymine in DNA; uracil in RNA
        m       a or c
        r       a or g
        w       a or t
        s       c or g
        y       c or t
        k       g or t
        v       a or c or g; not t
        h       a or c or t; not g
        d       a or g or t; not c
        b       c or g or t; not a
        n       a or c or g or t
