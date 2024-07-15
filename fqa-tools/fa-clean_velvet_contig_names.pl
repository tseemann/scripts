#!/usr/bin/env perl
#
#

use strict;
use Bio::SeqIO;

my $sii = Bio::SeqIO->new(-file => "$ARGV[0]", -format => 'Fasta');

my $sio = Bio::SeqIO->new(-file => ">$ARGV[1]", -format => 'Fasta');

while (my $s = $sii->next_seq()){
	my $id = $s->id;
	$id =~ s/_length.*//;
	$s->id($id);
	$sio->write_seq($s);
}
