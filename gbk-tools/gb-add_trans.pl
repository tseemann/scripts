#!/usr/bin/env perl

#  gb-ass_trans.pl
#  
#  Copyright 2013 Simon Gladman
#  Victorian Bioinformatics Consortium
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

use warnings;
use strict;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";
use Bio::SeqIO;
use Bio::Seq;

my(@Options, $verbose, $table);
setOptions();


#----------------------------------------------------------------------
# Reads a GenBank file, checks to see if each CDS has a translation
# annotation and if not adds the correct one.

my $in   = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'GenBank');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'GenBank');
while (my $seq = $in->next_seq) {
    for my $f($seq->all_SeqFeatures){
        if ($f->primary_tag eq "CDS"){
            my $s = $f->seq;
            my $p = $s->translate(-codontable_id=>$table);
            #print STDERR $p->seq . "\n";
            unless($f->has_tag("translation")){
                $f->add_tag_value("translation", $p->seq);
            }
        }
    }
  $out->write_seq($seq);
} 

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
	use Getopt::Long;

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
		{OPT=>"t|table=i", VAR=>\$table, DEFAULT=>11, DESC=>"The codon translation table to use."},
    );

	(@ARGV < 1) && (usage());

	&GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

	# Now setup default values.
	foreach (@Options) {
		if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
		${$_->{VAR}} = $_->{DEFAULT};
		}
	}
}

sub usage {
	print "Usage: $0 [options] gbk_file > new.gbk\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	exit(1);
}
 
#----------------------------------------------------------------------
