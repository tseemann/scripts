#!/usr/bin/env perl
#
#

use strict;
use Bio::SeqIO;
use Bio::Seq;

sub by_num { $b <=> $a }

my $scaffile = $ARGV[0];
my $reportfile = $ARGV[1];

my %scaffs = ();
my %reports = ();

open IN, $reportfile;


my $count = 0;

while(<IN>){
    if(/\w+\s+\d+/){
        my @tmp = split /\s+/, $_;
        $reports{$tmp[0]}->{$tmp[1]}->{'Type'} = $tmp[2];
        $reports{$tmp[0]}->{$tmp[1]}->{'Old'} = $tmp[3];
        $reports{$tmp[0]}->{$tmp[1]}->{'New'} = $tmp[4];
        $count ++;
    }
}

print STDOUT "Loaded $count edit requests from $reportfile\n";

my $scio = Bio::SeqIO->new(-file => $scaffile, -format => 'Fasta');

while(my $seq = $scio->next_seq()){
    $scaffs{$seq->id} = $seq->seq();
    print STDOUT "Loaded ". $seq->id() . "\n";
}

foreach my $key (sort keys %reports){
    print STDOUT "Working on $key..\n";
    #print STDOUT "Old seq: " . $scaffs{$key} . "\n";
    my %edits = %{$reports{$key}};
    foreach my $pos(sort by_num keys %edits){
        print STDOUT "Working on position $pos..\n";
        print STDOUT $edits{$pos}->{'Type'} . "\t" . $edits{$pos}->{'Old'} . "\t" . $edits{$pos}->{'New'} . "\n";
        my $t = $edits{$pos}->{'Type'};
        if($t =~ /subs/i){
            #substitution
            print STDOUT "Its a substitution!\n";
            substr($scaffs{$key}, ($pos-1), 1) = $edits{$pos}->{'New'};
        }
        elsif($t =~ /delet/i){
            #deletion!
            print STDOUT "Its a deletion!\n";
            substr($scaffs{$key}, ($pos-1), 1) = "";
        }
        elsif($t =~ /insertion/i){
            #insertion
            print STDOUT "Its an insertion!\n";
            substr($scaffs{$key}, ($pos-1), 0) = $edits{$pos}->{'New'};
        }
        #print STDOUT "New seq: " . $scaffs{$key} . "\n";
    }
}

#now output everything via a bio::seqio object...

my $outfile = $scaffile . "_shrimped.fna";

my $seqout = Bio::SeqIO->new(-file => ">$outfile", -format => 'Fasta');

foreach my $key (sort keys %scaffs){
    print STDOUT "Outputting new version of $key\n";
    my $seq = Bio::Seq->new(-id => $key, -seq => $scaffs{$key});

    $seqout->write_seq($seq);
}
