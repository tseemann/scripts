#!/usr/bin/env perl
use strict;
use warnings;

#    Copyright (C) 2011 Torsten Seemann <torsten@seemann.id.au>
#   
#	http://bioinformatics.net.au
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

my $FASTQ_OFFSET = 33;    # 33=sanger 64=illumina(<1.8)
my $AFG_OFFSET = 60;

my $iid=0;
while (my $eid = <ARGV>) {

  die "bad fastq ID line '$eid'" unless $eid =~ m/^\@/;
  $eid = substr $eid, 1;  # remove '@'

  my $seq = scalar(<ARGV>);
  chomp $seq;

  my $id2 = scalar(<ARGV>);
  die "bad fastq ID2 line '$id2" unless $id2 =~ m/^\+/;

  my $qlt = scalar(<ARGV>);
  chomp $qlt;

  $iid++;
  print "{RED\n";
  print "iid:$iid\n";
  print "eid:$eid";  # already has \n
  print "seq:\n";  
  for (my $i=0; $i < length($seq); $i+=60) {
    print substr($seq, $i, 60), "\n";
  }
  print ".\n";
  print "qlt:\n";  
  my @q = split m//, $qlt;
  @q = map { chr( ord($_)-$FASTQ_OFFSET + $AFG_OFFSET ) } @q;
  for (my $i=0; $i < @q; $i++) {
    print $q[$i];
    print "\n" if ($i+1) % 60 == 0;
  }
  print ".\n";
  print "}\n";
}


#{RED
#iid:1
#eid:1
#seq:
#BCAAAAAGAACTGACTGCCCACCAACTAGAATTCGTATTGGGTATAAAACCTATTACGAAT
#TAATGCAGAATCCT
#.
#qlt:
#CAAAAAGAACTGACTGCCCACCAACTAGAATTCGTATTGGGTATAAAACCTATTACGAAT
#TAATGCAGAATCCT
#.
#}
