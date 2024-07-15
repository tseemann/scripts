#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use File::Temp;
use Bio::SearchIO;
use Data::Dumper;

my(@Options, $verbose, $informat, $outformat, $hmmerdb, $blastdb, $evalue, $cores);
setOptions();

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>$informat);
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>$outformat);

while (my $seq = $in->next_seq) 
{
  for my $f ($seq->get_SeqFeatures) 
  {
    next unless $f->primary_tag eq 'CDS';
    my $id = ($f->get_tag_values('locus_tag'))[0];
    next unless $id;
    my $desc = ($f->get_tag_values('product'))[0] || '';
    next unless !$desc or $desc =~ m/hypothetical/;
    my $prod = identify_protein_HMMER($f->seq->translate);
    print STDERR "$id | $desc => ", ($prod || '[none found]'), "\n";
    if ($prod) {
      $f->remove_tag('product');
      $f->add_tag_value('product', $prod);
    }
  }
  $out->write_seq($seq);
}
print STDERR "\nDone\n";

#----------------------------------------------------------------------

sub identify_protein_HMMER {
  my($pep) = @_;
  my $fh = File::Temp->new();
  my $out = Bio::SeqIO->new(-file=>">$fh", -format=>'fasta');
  $out->write_seq($pep);
  my $cmd = "hmmscan --cpu $cores -E $evalue $hmmerdb $fh";
  print STDERR "Running: $cmd\n";
  my $report = Bio::SearchIO->new(-file=>"$cmd |", -format=>'hmmer3') or return;
  while (my $result = $report->next_result) {
    while (my $hit = $result->next_hit) {
      if ($hit->description !~ m/hypothetical|DUF/) {
        return  $hit->description;
      }
    }
  }
  return '';
}

#----------------------------------------------------------------------

sub identify_protein_BLAST {
  my($pep) = @_;
#  print Dumper($pep); exit;
  my $fh = File::Temp->new();
  my $out = Bio::SeqIO->new(-file=>">$fh", -format=>'fasta');
  $out->write_seq($pep);
  my $cmd = "blastp -query $fh -evalue $evalue -db $blastdb -num_threads $cores";
  print STDERR "Running: $cmd\n";
  my $report = Bio::SearchIO->new(-file=>"$cmd |", -format=>'blast') or return;
  while (my $result = $report->next_result) {
    while (my $hit = $result->next_hit) {
      if ($hit->description !~ m/hypothetical|DUF/) {
        return $hit->description;
      }
    }
  }
  return;
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose+",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"if=s",  VAR=>\$informat, DEFAULT=>'genbank', DESC=>"Input format: genbank embl"},
    {OPT=>"of=s",  VAR=>\$outformat, DEFAULT=>'genbank', DESC=>"Output format: genbank embl"},
    {OPT=>"blastdb=s",  VAR=>\$blastdb, DEFAULT=>'/bio/db/blast/refseq_microbial_protein', DESC=>"BLAST protein database"},
    {OPT=>"hmmerdb=s",  VAR=>\$hmmerdb, DEFAULT=>'/bio/db/hmmer3/PRK.hmm', DESC=>"BLAST protein database"},
    {OPT=>"evalue=f",  VAR=>\$evalue, DEFAULT=>1E-6, DESC=>"Similarity e-value cut-off"},
    {OPT=>"cores=i",  VAR=>\$cores, DEFAULT=>8, DESC=>"Threads to use"},
  );

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] old.gbk > new.gbk\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
