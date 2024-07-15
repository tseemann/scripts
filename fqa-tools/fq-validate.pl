#!/usr/bin/env perl

while ( not eof() ) {

  my $id = scalar(<>);
  chomp $id;
  die "bad id: $id" unless substr($id,0,1) eq '@';

  my $seq = scalar(<>);
  chomp $seq;
  die "bad seq: $seq" unless $seq =~ m/^[actgun]+$/i;

  my $div = scalar(<>);
  chomp $div;
  die "bad div: $div" unless substr($div,0,1) eq '+';
  
  die "id not match div:\nid=$id\ndv=$div" if length($div)>1 and substr($id,1) ne substr($div,1);

  my $qual = scalar(<>);
  chomp $qual;
  die "qual not match seq:\nid==$id\nseq=$seq\nqul=$qual" if length($seq) != length($qual);
  
}