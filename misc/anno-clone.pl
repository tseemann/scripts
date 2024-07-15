#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Temp qw(tempdir);
use File::Spec;

#    anno-clone.pl
#
#    Takes unannotated DNA sequence(s) + a related annotated DNA sequence
#    and transfers the annotation to the unannotated sequence(s).
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

my(@Options, $verbose, $protein, $informat, $outformat, $join, $evalue);
setOptions();
sub inform;

#$protein && problem("--protein not implemented yet");
$join && problem("--join not implemented yet");

require_exe('blastall');
require_exe('formatdb');

my $tmpdir = tempdir(CLEANUP => 1);

my $featname = File::Spec->catfile($tmpdir, 'features.fasta');
my $out = Bio::SeqIO->new(-file=>">$featname", -format=>'Fasta');
inform "Writing: $featname (fasta)";

my $refname = $ARGV[1];
inform "Reading: $refname ($informat)";
my $ref = Bio::SeqIO->new(-file=>$refname, -format=>$informat);
my %feat;
my $numfeat=0;
while (my $seq = $ref->next_seq) {
  for my $f ($seq->get_SeqFeatures) {
#    next if $f->primary_tag =~ m/^(source|gene)$/;
    next if $f->primary_tag eq 'source';
    $numfeat++;
    $feat{$numfeat} = $f;
    my $s = $f->spliced_seq;
    $s->display_name($numfeat);
    $s->desc('');
    $out->write_seq($s);
  }
}

my $refdbname = File::Spec->catfile($tmpdir, 'refdb');
runcmd("formatdb -i $ARGV[0] -n $refdbname -p F -l /dev/null");
#system("ls -l $tmpdir");

my $blsname = "$featname.blast";
runcmd("blastall -p blastn -i $featname -o $blsname -d $refdbname -F F -v 1 -b 1 -e $evalue");
#system("ls -l $tmpdir");
#system("cat $blsname");

my $ctgname = $ARGV[0];
inform "Loading: $ctgname";
my $ctg = Bio::SeqIO->new(-file=>$ctgname, -format=>'Fasta');
my %ctg;
while (my $seq = $ctg->next_seq) {
  $ctg{$seq->id} = $seq;
}

inform "Reading: $blsname";
my $bls = Bio::SearchIO->new(-file=>$blsname, -format=>'blast');
while (my $res = $bls->next_result) {
#  inform "BLAST: ",$res->query_name;
  while (my $hit = $res->next_hit) {
    while (my $hsp = $hit->next_hsp) {
      my $f = Bio::SeqFeature::Generic->new(-gff_string=>$feat{$res->query_name}->gff_string);
      $f->seq_id( $hit->name);
      $f->start( $hsp->start('hit') );
      $f->end( $hsp->end('hit') );
      $f->strand( $hsp->strand('hit') );
#      print $f->gff_string,"\n";
      $ctg{$hit->name}->add_SeqFeature($f);
    }
  }
}

inform "Writing: STDOUT (genbank)";
my $anno = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>$outformat);
for my $ctg (values %ctg) {
  $anno->write_seq($ctg);
}


#-------------------------------------------------------------------------

sub runcmd {
  inform "Running: @_";
  system(@_)==0 or problem("system(@_) failed: $?");
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

sub require_exe {
  my($bin) = shift;
  for my $dir (File::Spec->path) {
    my $exe = File::Spec->catfile($dir, $bin);
    return $exe if -x $exe; 
  }
  return;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"p|protein!",  VAR=>\$protein, DEFAULT=>0, DESC=>"Use protein matching, not DNA"},
    {OPT=>"if|informat=s",  VAR=>\$informat, DEFAULT=>'genbank', DESC=>"Annotated reference input format"},
    {OPT=>"of|outformat=s",  VAR=>\$outformat, DEFAULT=>'genbank', DESC=>"Annotated output format"},
    {OPT=>"j|join!",  VAR=>\$join, DEFAULT=>0, DESC=>"Output single pseudo-molecule"},
    {OPT=>"e|evalue=f",  VAR=>\$evalue, DEFAULT=>1E-3, DESC=>"BLAST e-value cutoff"},
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] <contigs.fasta> <reference.gbk>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
