#!/usr/bin/env perl

#I don't know if there's a way of doing this in Artemis (if not, it would
#certainly be a useful feature to add), but here is a simple script to do
#what you want. You'll need to have Perl and BioPerl installed. The usage is
#(assuming you save the program as an executable file called "nucrestart",
#and your sequence is in the file "seq.embl"): nucrestart seq.embl 1000
#--format embl > newseq.embl The file "newseq.embl" should now contain a
#version of the sequence with the old base 1000 as the new base 1. You can
#also include the option "--revcom" to reverse complement the sequence. All
#the coordinates of the features on the sequence are adjusted appropriately.

use warnings FATAL=>qw(all);
use strict;
use Bio::SeqIO;
use Bio::SeqUtils qw(cat trunc_with_features revcom_with_features);
use Getopt::Long;

my $revcom=0;
my $format='genbank';
GetOptions('revcom!'=>\$revcom, 'format=s'=>\$format) or die 'Died';

die "Changes the start position of the sequence (assumes circularity).\nUsage: $0 seq_file 100\n" unless defined $ARGV[1];

my $seq=Bio::SeqIO->new(-format=>$format, -file=>$ARGV[0])->next_seq;
die "Error - New start not numeric\n" if $ARGV[1]=~/\D/;
die "Error - New start not in sequence\n" if $ARGV[1]>$seq->length or $ARGV[1]==0;
my $first=Bio::SeqUtils->trunc_with_features($seq, $ARGV[1], $seq->length);
my $second=Bio::SeqUtils->trunc_with_features($seq, 1, $ARGV[1]-1);
Bio::SeqUtils->cat($first, $second);
$first=Bio::SeqUtils->revcom_with_features($first) if $revcom;
my $out=Bio::SeqIO->newFh(-format=>$format);
print $out $first;
