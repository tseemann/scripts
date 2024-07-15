#!/usr/bin/env perl

use strict;

die "program csfasta shift\nif shift is 1, the first base is omitted in the output\n" if (@ARGV < 1);

my $csfasta_file = $ARGV[0];
my %colourspace = colour_space();
my $shift = 0;
my %sequences = read_csfasta( $csfasta_file );
$shift = $ARGV[1] if( @ARGV == 2 );

foreach my $key ( keys %sequences ){
    print ">$key\n"; 
    my $seq = $sequences{ $key };
    my @letters = split( //, $seq );
    my $first_base = $letters[0];
    
    for( my $i = 1; $i < @letters ; $i++ ){
	
	my $colour = $letters[$i];
	my $encoding = $first_base.$colour;
	$first_base = $colourspace{ $encoding };
	$letters[ $i ] = $first_base;    
    }
    
    shift( @letters ) if( $shift );
    $" = "";
    print "@letters\n";

}


sub colour_space{

    my %hash = (
	"A0" => "A",
	"C0" => "C",
	"G0" => "G",
	"T0" => "T",
	"A1" => "C",
	"C1" => "A",
	"G2" => "A",
	"A2" => "G",
	"A3" => "T",
	"T3" => "A",
	"C2" => "T",
	"T2" => "C",
	"C3" => "G",
	"G3" => "C",
	"G1" => "T",
	"T1" => "G" );
    return %hash;


}




sub read_csfasta{
    my ( $csfasta_file) = @_;
    my %sequences;
    my $first = 1;
    
    my $header;
    my $sequence;
    
    open(CSFASTA, $csfasta_file);
    while( my $line = <CSFASTA> ){
        
	chomp $line;    
	
	if( ! $first ){     # $sequence and header are not initialized in first iteration 
	    
	    if(  $line =~ /\>/   ){         # encounters a new sequence             
		
		$sequences{$header} = $sequence;
		$header = $line; $header =~ tr/\>//d;  # read a new header line
		$sequence = "";
	    }
	    
	    else{ $sequence .= $line;       # concatenate sequences
	    }
	}
	else{                                       # read the header line
	    if(  $line =~ /\>/   ){ 
	    
		$header = $line; $header =~ tr/\>//d;$sequence = "";
		
	    }
	}       
	$first = 0;
    }
    close( CSFASTA );
    
    ## the last one
    $sequences{$header} = $sequence;    
    
    return %sequences;
}
