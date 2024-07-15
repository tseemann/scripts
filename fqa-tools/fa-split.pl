#!/usr/bin/env perl

use strict;
use Bio::SeqIO;

my(@Options, $verbose, $template, $num, $trunc, $raw);
setOptions();

my $in  = Bio::SeqIO->new(-fh => \*ARGV, -format => 'Fasta');
my $fcount = 0;
my $scount = $num+1;
my $out;
my $outname;

if ($template =~ m{/}) {
  die "Slash '/' not allowed in --template" if $template =~ m{/};
}
if ($template =~ m/%i/ and $num != 1) {
  die "Must set --num=1 if using %i in --template";
}

my @fofn;

while (my $seq = $in->next_seq) {
  if ($scount >= $num) {
    $fcount++;
    $scount=0;
    $outname = $template;
    my $id = $seq->display_id;
    $outname =~ s/%i/$id/g;
    $outname =~ s/%d/$fcount/g;
    $outname =~ s/[^\w\-\.]/_/g;
    $outname = sprintf $outname, $fcount;
    die "Output '$outname' has slash '/' in it!" if $outname =~ m{/};
    push @fofn, $outname if $raw;
    $out = Bio::SeqIO->new(-file => ">$outname", '-format' => ($raw ? 'raw' : 'Fasta') )
      or die "coult not write to '$outname'";
  }
  $seq = $seq->trunc(1,$trunc) if $trunc and $seq->length > $trunc;
  $out->write_seq($seq);
  $scount++;
  if ($verbose) {
    print STDERR "[$fcount.$scount] Wrote '",$seq->id,"' to '$outname'.\n";
  }
}

if ($raw) {
  print STDERR "Writing 'fofn' file of filenames for GAP4\n";
  open FOFN, '>', 'fofn';
  print FOFN (map { "$_\n" } @fofn);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!", VAR=>\$verbose, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"template=s", VAR=>\$template, DEFAULT=>'%i.seq', 
    	DESC=>'Pattern to name files: %d = file_no, %i = fasta_id'},
    {OPT=>"n|num=i", VAR=>\$num, DEFAULT=>1, DESC=>"Seqs per file"},
    {OPT=>"raw!", VAR=>\$raw, DEFAULT=>0, DESC=>"Raw output with fofn for GAP4"},
    {OPT=>"trunc=i", VAR=>\$trunc, DEFAULT=>undef, 
    	DESC=>"Truncate each sequence to this many symbols"}
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
  print "Usage: $0 [options] file1.fna [ file2.fna ... ]\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
