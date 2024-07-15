#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $debug, $format, $idtag, $prot, $ftype, $desctag, $species);
setOptions();

#print STDERR "Opening STDOUT ...\n";
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

#print STDERR "Loading STDIN ...\n";
my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>$format);

while (my $seq = $in->next_seq) {
  print STDERR "\rParsing: ",$seq->display_id;
  my $genspec = $species ? ' ['.($seq->species->binomial('FULL') || 'Genus species').']'
                         : '';
  my $counter=0;
  for my $f ($seq->get_SeqFeatures) {
#    print STDERR "\r", $seq->display_id, " ", $f->primary_tag, " ", $f->location->to_FTstring,"\n" if $debug;
    next if $f->primary_tag eq 'source';
    next if $ftype and $f->primary_tag ne $ftype;
    $counter++;
    print STDERR "\rProcessing: ", $seq->display_id, " | $counter";
    my $id = $f->has_tag($idtag) ? ($f->get_tag_values($idtag))[0] : "SEQ$counter";
    my $desc = $f->has_tag($desctag) ? ($f->get_tag_values($desctag))[0] : "not_available";
    $desc .= $genspec;
#    my $s = $seq->trunc( $f->location );
    my $s = $f->spliced_seq;  # don't forget eukaryotes!
    # HANDLE CODON START FOR FUZZY FEATURES!
    if ($f->has_tag('codon_start')) {
      my($cs) = $f->get_tag_values('codon_start');
      if ($cs != 1) {
        print STDERR "/codon_start=$cs - trimming mRNA!";
        $s = $s->trunc($cs, $s->length);
      }
    }
    #END
    $s = $s->translate if $prot;
    $out->write_seq(Bio::Seq->new(-display_id=>$id, -seq=>$s->seq, -desc=>$desc));
  }
}
print STDERR "\nDone\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"format=s",  VAR=>\$format, DEFAULT=>'genbank', DESC=>"Input format"},
    {OPT=>"ftype=s",  VAR=>\$ftype, DEFAULT=>'CDS', DESC=>"Which feature type"},
    {OPT=>"idtag=s",  VAR=>\$idtag, DEFAULT=>'locus_tag', DESC=>"What tag to use as Fasta ID"},
    {OPT=>"desctag=s",  VAR=>\$desctag, DEFAULT=>'product', DESC=>"What tag to use as Fasta DESC"},
    {OPT=>"species!",  VAR=>\$species, DEFAULT=>0, DESC=>"Include [Genus species] in desc"},
#    {OPT=>"trunc=i",  VAR=>\$trunc, DEFAULT=>0, DESC=>"+ve N-term, -ve C-term"},
    {OPT=>"prot!",  VAR=>\$prot, DEFAULT=>0, DESC=>"As protein?"},
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
  print "Usage: $0 [options] < file.gbk > file.fna\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
