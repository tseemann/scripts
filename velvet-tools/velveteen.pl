#!/usr/bin/perl -w
use strict;
use Fatal;
use File::Spec;
use Data::Dumper;
use Time::Piece;

#    velveteen.pl
#
#    Provides a DWIM ("do what I mean") interface to velveth/velvetg
#
#    Copyright (C) 2009 Torsten Seemann
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

if (@ARGV < 1 or $ARGV[0] =~ m/^--?h/i) {
  print<<"EOF";
NAME
  velveteen - makes it easy to run velveth/velvetg for simple cases
            - the equivalent of "runAssembly *.sff *.fasta" for 454 data
SYNOPSIS
  $0 <k-value | single.fa[.gz] | pair.fa[.gz] | single.fq[.gz] | pair.fq[.gz]>
EXAMPLES
  $0 reads.fasta
  $0 s_1_sequence.txt 29 454Contigs.fna
  $0 31 *.fq 
DESCRIPTION
  . parameters can be in any order
  . auto-detects if velveth and velvetg are installed
  . auto-detects fasta or fastq files
  . assumes gzipped files if name ends in .gz
  . auto-detects single or paired reads (shuffled with /1 /2 or _R1 _R2 notation)
  . auto-chooses paired-mode if read file name contains word "pair" in it.
  . defaults to k=31 unless an integer specified
  . output folder is ./velvet-KK-YYYY-MM-DD-HH-MM-DD/
  . uses -exp_cov auto and -cov_cutoff auto and auto insert length estimation
BUGS
  . does not manage more than one category/fileset correctly
  . assumes sequences are not multi-lines in the .fq/.fa files
  . output folder not mutable
AUTHOR
  Torsten Seemann
URL
  http://www.bioinformatics.net.au/
EOF
  exit 1;
};

for my $exe ('velveth', 'velvetg', 'zcat') {
  my $path = require_exe($exe) or problem("$exe not found in PATH");
  inform("$exe ... using $path");
}
my $CATS = get_setting("velveth |", 'CATEGORIES\s*=\s*(\d+)');
inform("categories = $CATS");
my $MAXK = get_setting("velveth |", 'MAXKMERLENGTH\s*=\s*(\d+)');
inform("maxkmerlength = $MAXK");

my $k=31;
#inform("default k = $k");
my @file;

for my $v (@ARGV) {
  if ($v =~ /^\d+$/ and $v <= $MAXK and $v >= 2) {
    $k = $v;
    inform("set k = $k");
  }
  elsif (-r $v) {
    push @file, $v;
    inform("added input file = $v");
  }
  else {
    inform("ignoring parameter '$v'");
  }
}

problem("No readable sequence input files provided") unless @file;

my $t = localtime;
my $dir = "./velvet-$k-".$t->ymd('.').'-'.$t->hms('.');
inform("output directory = $dir");

my $vhcmd = "velveth $dir $k -create_binary";
my %file;
for my $file (@file) {
  push @{ $file{ guess_input($file) } }, $file;
}
for my $type (sort keys %file) {
  inform("$type = @{$file{$type}}");
  $vhcmd .= " $type @{$file{$type}}";
}

#  -clean <yes|no> : remove all the intermediary files which are useless for recalculation (default : no)
#  -very_clean <yes|no> : remove all the intermediary files (no recalculation possible) (default: no)
    

runcmd($vhcmd);
runcmd("velvetg $dir -clean yes -exp_cov auto -cov_cutoff auto 1> $dir/velvetg.stdout 2> $dir/velvetg.stderr");
inform("results in $dir");

#-------------------------------------------------------------------------

sub guess_input {
  my($fname) = @_;
  my $gzipped = $fname =~ m/gz$/i ? '.gz' : '';
  my $cat = '';
  open FH, "zcat -f $fname |";
  my @line = map { scalar(<FH>) } (1..8); # read first 8 lines
  chomp @line;
  for my $i (0..7) {
#    print STDERR "[$i] $line[$i]\n";
  }
  # Old Illumina reads had /1 and /2 for the F and R read IDs suffix
  # But OLD >1.8 have the SAME ID, and put 1: 2: in the description instead!
  
  if ($line[0] =~ m/^\s*>/) {
    if ( $fname=~m/pair/i or ($line[0] =~ m{/1\s*$} && $line[2] =~ m{/2\s*$}) or ($line[0] eq $line[2]) ) {
      return "-shortPaired -fasta$gzipped";
    }
    return "-short -fasta$gzipped";
  }
  elsif ($line[0] =~ m/^\@/) {
    if ( $fname=~m/pair/i or ($line[0] =~ m{/1\s*$} and $line[4] =~ m{/2\s*$}) or ($line[0] eq $line[4]) ) {
      return "-shortPaired -fastq$gzipped";
    }
    return "-short -fastq$gzipped";
  }
}

#-------------------------------------------------------------------------

sub get_setting {
  my($fio,$regexp) = @_;
  open FH, $fio;
  while (<FH>) {
    return $1 if m/$regexp/;
  }
  return;
}


#-------------------------------------------------------------------------

sub inform {
  my @s = map { defined $_ ? $_ : '(undef)' } @_;
  my $t = localtime;
  print STDERR "[$t] ",@s,"\n";
}

#-------------------------------------------------------------------------

sub problem {
  inform 'ERROR: ', @_;
  exit -1;
}

#-------------------------------------------------------------------------

sub runcmd {
  inform "Running: @_";
  system(@_)==0 or problem("system(@_) failed: $?");
}

#-------------------------------------------------------------------------

sub require_exe {
  my($bin) = shift;
  for my $dir (File::Spec->path) {
    my $exe = File::Spec->catfile($dir, $bin);
    return $exe if -x $exe; 
  }
  return;
}

#----------------------------------------------------------------------
