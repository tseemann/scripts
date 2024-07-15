#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;
use Text::CSV;

my(@Options, $verbose, $format);
setOptions();

my $csv = Text::CSV->new;

open my $fh, "zcat -f @ARGV |";
my $gbk = Bio::SeqIO->new(-fh=>$fh, -format=>$format);

print "SOURCE,GENE,LOCUS_TAG,PRODUCT,CDNA,TRANSLATION,GO_TERMS\n";

while (my $seq = $gbk->next_seq) {
  print STDERR "Processing: ",$seq->id,"\n";
  my @row = ($seq->id);
  for my $f ($seq->get_all_SeqFeatures) {
#    print STDERR "\tFeature: ",$f->primary_tag,"\n";
    if ($f->primary_tag eq 'gene') {
      push @row, tag($f,'gene','');
      push @row, tag($f,'locus_tag','');
      push @row, tag($f,'note','hypothetical');
    }
    elsif ($f->primary_tag eq 'mRNA') {
      push @row, $f->spliced_seq->seq;
    }
    elsif ($f->primary_tag eq 'CDS') {
      push @row, tag($f,'translation','');
      push @row, tag($f,'db_xref','', '^GO:');
      $csv->combine(@row) || die $csv->error_diag;
      print $csv->string,"\n";
      @row=($seq->id);
    }
  }
}

sub tag {
  my($f, $tag, $def, $patt) = @_;
  if ($f->has_tag($tag)) {
    my(@val) = $f->get_tag_values($tag);
    @val = grep(m/$patt/ , @val) if $patt;
    return join(',',@val);
  }
  return $def;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"f|format=s",  VAR=>\$format, DEFAULT=>'genbank', DESC=>"genbank | embl"},
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
  print "Usage: $0 [options] file1.gbk.gz [file2.gbk ...] > output.csv\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__
     gene            14403..20106
                     /gene=ENSPVAG00000016071
     mRNA            join(14403..14786,16351..16559,17637..17697,17828..17923,
                     18378..18521,18523..18543,18801..18926,19200..19420,
                     20076..20106)
                     /gene="ENSPVAG00000016071"
                     /note="transcript_id=ENSPVAT00000016071"
     CDS             join(14403..14786,16351..16559,17637..17697,17828..17923,
                     18378..18521,18523..18543,18801..18926,19200..19420,
                     20076..20106)
                     /gene="ENSPVAG00000016071"
                     /protein_id="ENSPVAP00000015170"
                     /note="transcript_id=ENSPVAT00000016071"
                     /db_xref="Ens_Hs_transcript:ENST00000293670"
                     /db_xref="Ens_Hs_translation:ENSP00000293670"
                     /db_xref="GO:GO:0004867"
                     /db_xref="GO:GO:0005198"
                     /db_xref="GO:GO:0045095"
                     /translation="MTCGFSNVGCGFGPRNFSCASACGPRPGRCCITAAPYRGISCYR
                     GLTGGFGSRSVCGGFRSGSCGRSFGYRSGGVCGPSPPCITTVSVNESLLTPLNLEIDP
                     NAQCVKHEEKEQIKCLNSRFAAFIDKVRFLEQQNKLLETKWQFYQNRKCCESNLEPLF
                     EGYIETLRREAECVEADSGRLASELNHVQEVLEGYKKKYEEEVALRATAENQFVALKK
                     DVDCAYLRKSDLEANVEALTQEIDFLRRLYEEXXXXXXXXXXXXXXXXXXXXXXXXXX
                     XXXXXXXXXXXXXXXXXXXXXXESWYRSKCEEMKATVIRHGETLRRTKEEINELNRMI
                     QRLTAEVENAKCQNTKLEAAVTQSEQQGETAISDARCKLAELEAALQKAKQDMACLLR
                     EYQEVMNSKLGLDIEIATYRRLLEGEEQRLCEGIGAVNV"
