#!/usr/bin/env perl

#    Copyright (C) 2011 Torsten Seemann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


use strict;
use Bio::SeqIO;

my(@Options, $verbose, $ref_fn, $primers_fn, $maxsize, $exclude);
setOptions();

-r $ref_fn or die "please specify -r reference file";
my $ref_io = Bio::SeqIO->new(-file=>$ref_fn, -format=>'Fasta');

-r $primers_fn or die "please specify -p primers file";

my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

while (my $ref = $ref_io->next_seq) 
{
  print STDERR "Ref: ", $ref->display_id,"\n";
  my @pair;
  my @name;
  my $primers_io = Bio::SeqIO->new(-file=>$primers_fn, -format=>'Fasta');
  
  while (my $p = $primers_io->next_seq) 
  {
    # check primer
    print STDERR "Checking ",$p->display_id, " fwd: ", $p->seq, " - ";
    my $fp = find_it($p->seq, $ref->seq);
    print STDERR ($fp ? "$fp\n" : "no\n");
    
    # check revcom primer
    print STDERR "Checking ",$p->display_id, " rev: ", $p->revcom->seq, " - ";
    my $rp = find_it($p->revcom->seq, $ref->seq);
    $rp += $p->length-1 if $rp; # correct position for fwd sense if found
    print STDERR ($rp ? "$rp\n" : "no\n");

    # keep track of positions found    
    push @pair, ($fp || $rp);
    push @name, $p->display_id;

    # we have now processed two primers so let's make sense of them    
    if (@pair == 2) 
    {
      my $product = '';
      if (abs($pair[0] - $pair[1]) > $maxsize) {
        # hmmm, too big, bad primers?
        print STDERR "\n**WARNING** product too large for: @name | @pair\n\n";
      }
      elsif ($pair[0] < $pair[1]) {
        # fwd product
        $product = $ref->subseq($pair[0], $pair[1]);
      }
      else {
        # rev product
        $product = $ref->trunc($pair[1]-60, $pair[0]+60)->revcom->seq;
      }
      
      # output product if we have one
      if ($product) {
	my $desc = length($product)."bp $pair[0]..$pair[1] ".$ref->display_id;
        print STDERR "Found product: $desc\n";
	$out->write_seq( 
          Bio::Seq->new(
	    -id   => join('+', @name),
	    -desc => $desc,
	    -seq  => $product,
	  )
	);
      }
      
      # reset ready for next pair
      @pair = ();
      @name = ();
    }
  }

}
print STDERR "Done.\n";


#----------------------------------------------------------------------

sub find_it {
  my($p, $r) = @_;
  my $pos = index($r, $p);
  return $pos+1 if $pos >= 0;
  my $patt = string_ed1($p);
  $r =~ m/$patt/ig;
  $pos = pos($r) || -1;
  return $pos+1;
}

#----------------------------------------------------------------------

sub string_ed1 {
  my($s) = @_;
  my $L = length $s;
  my @patt = $s;
  for my $i (0 .. $L-1) {
    my $p = $s;
    substr $p, $i, 1, '.';
    push @patt, $p;
  }
  return '('.join('|',@patt).')';
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"r|reference=s",  VAR=>\$ref_fn, DEFAULT=>'', DESC=>"Reference sequence, fasta"},
    {OPT=>"p|primers=s",  VAR=>\$primers_fn, DEFAULT=>'', DESC=>"Primer pairs, fasta, layout L1 R1 L2 R2 ..."},
    {OPT=>"m|maxsize=i",  VAR=>\$maxsize, DEFAULT=>10_000, DESC=>"Maximum expected product size"},
    {OPT=>"x|exclude!",  VAR=>\$exclude, DEFAULT=>0, DESC=>"Exclude the primer sequences"},
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
  print "Usage: $0 [options] -r ref.fa -p primers.fa > products.fa\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__

Hi Torst

that is very nice of you.  I was just in the process of working out BLAST and bioperl to pull out that region from the fasta.  The application I have in mind is in silico MLST, so lets start with the PA01 genome

http://www.ncbi.nlm.nih.gov/nucleotide/NC_002516?&report=fasta

and the MLST primers

>acsA-F
ACCTGGTGTACGCCTCGCTGAC
>acsA-R
GACATAGATGCCCTGCCCCTTGAT
>aroE-F
TGGGGCTATGACTGGAAACC
>aroE-R
TAACCCGGTTTTGTGATTCCTACA
>guaA-F
CGGCCTCGACGTGTGGATGA
>guaA-R
GAACGCCTGGCTGGTCTTGTGGTA
>mutL-F
CCAGATCGCCGCCGGTGAGGTG
>mutL-R
CAGGGTGCCATAGAGGAAGTC
>nuoD-F
ACCGCCACCCGTACTG
>nuoD-R
TCTCGCCCATCTTGACCA
>ppsA-F
GGTCGCTCGGTCAAGGTAGTGG
>ppsA-R
GGGTTCTCTTCTTCCGGCTCGTAG
>trpE-F
GCGGCCCAGGGTCGTGAG
>trpE-R
CCCGGCGCTTGTTGATGGTT

cheers
jason
