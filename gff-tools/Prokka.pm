package Prokka;

use Exporter qw(import);
@EXPORT_OK = qw(make_id require_exe inform problem num_cpu read_gff3 write_gff3 sort_gff3);

use constant GFF_VERSION => 3;

use strict;
use warnings;
use Time::Piece;
use File::Spec;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use Carp;

#-------------------------------------------------------------------------

sub make_id {
  my($prefix, $value, $width) = @_;
  $width ||= 5;
  sprintf "%s%0*d", $prefix, $width, $value;
}

#-------------------------------------------------------------------------

sub inform {
  my @s = map { defined $_ ? $_ : '(undef)' } @_;
  my $t = localtime;
  print STDERR "[$t] ",@s,"\n";
}

#-------------------------------------------------------------------------

sub problem {
  inform @_;
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

#-------------------------------------------------------------------------

sub num_cpu {
  if ($^O =~ m/linux/i) {
    my($num) = qx(grep -c ^processor /proc/cpuinfo);
    chomp $num;
    return $num if $num =~ m/^\d+/;
  }
  return 1;
}

#-------------------------------------------------------------------------

sub BioSeqFeat_cmp {
  $a->seq_id cmp $b->seq_id 
  or 
  $a->start <=> $b->start
  or
  $a->length <=> $b->length
}

#-------------------------------------------------------------------------

sub sort_gff3 {
  my($gff) = @_;
  $gff->{'features'} = [ sort BioSeqFeat_cmp @{$gff->{'features'}} ];
}

#-------------------------------------------------------------------------

sub read_gff3 {
  my(@filename) = @_;
  my %gff;
  @filename = ('-') unless @filename;

  for my $gff_file (@filename) {
    inform "Loading: $gff_file";
    my $fio = Bio::Tools::GFF->new(-file=>$gff_file, -gff_version=>GFF_VERSION);
    $fio->ignore_sequence(0); 
    $fio->features_attached_to_seqs(1);
    while (my $f = $fio->next_feature) {
      push @{$gff{'features'}}, $f;
    }
#    $gff{'features'} = [ sort BioSeqFeat_cmp @{$gff{'features'}} ];
#    inform Dumper \%gff; exit;
    for my $seq ($fio->get_seqs) {
      $gff{'sequences'}{$seq->display_id} = $seq;
    }
  }
  return \%gff;
}

#-------------------------------------------------------------------------

sub write_gff3 {
  my($fh, $gff) = @_;

  my $gff_factory = Bio::Tools::GFF->new(-gff_version=>GFF_VERSION);

  print $fh "##gff-version ", GFF_VERSION, "\n";
  
  for my $id (keys %{$gff->{'sequences'}}) {
    print $fh "##sequence-region $id 1 ", $gff->{'sequences'}{$id}->length, "\n";
  }
  
#  for my $f (sort BioSeqFeat_cmp @{$gff->{'features'}}) {
  for my $f (@{$gff->{'features'}} ) {
    print $fh $f->gff_string($gff_factory),"\n";
  }

  if (scalar keys %{$gff->{'sequences'}}) {
    print $fh "##FASTA\n";
    my $seqio = Bio::SeqIO->new(-fh=>$fh, -format=>'fasta');
    $seqio->write_seq( sort { $a->display_id cmp $b->display_id } values %{$gff->{'sequences'}} );
  }
}

#-------------------------------------------------------------------------

1;
