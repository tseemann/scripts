#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Bio::SeqIO;

my(@Options, $debug, $format, $key, $ignore_csv, $keep_csv, $noheader, 
             $empty, $boolean, $sep, $revsuffix, $ftype);
setOptions();

my %ignore = (map { ($_ => 1) } split m/,/, $ignore_csv);
delete $ignore{$key};

my %keep = (map { ($_ => 1) } split m/,/, $keep_csv);
$keep{$key} = 1;

my %feat;
my %header = ($key => 9E9);
my $count=0;

my $in = Bio::SeqIO->new(-fh=>\*ARGV, -format=>$format);
while (my $seq = $in->next_seq) {
  print STDERR "Parsing: ",$seq->display_id,"\n";
  my $counter=0;
  for my $f ($seq->get_SeqFeatures) {
    $count++;
#    print STDERR "\rProcessing: ", $seq->display_id, " | $counter";
    next unless $f->has_tag($key);
    next unless $f->primary_tag =~ /$ftype/;
    my $nv = {
      start => $f->start,
      end => $f->end,
      strand => $f->strand,
      ftype => $f->primary_tag,
    };
    for my $tag ($f->get_all_tags) {
      next if $ignore{$tag};
      next unless $keep{$tag};
#      next if $tag eq $key;
      my($value) = $f->get_tag_values($tag);
      $value = $boolean if !defined($value) or $value eq '_no_value';
#      $value =~ s/$sep/_/g;  # ensure separator character not in a column value
      $nv->{$tag} = $value;
      $header{$tag}++;
    }  
    $feat{ $nv->{$key} } = $nv;
    if ($revsuffix) {
      my %dupe = %$nv;
      $dupe{$key} = $nv->{$key}.$revsuffix;
      $feat{ $nv->{$key}.$revsuffix } = \%dupe;
    }
  }
}
print STDERR Dumper(\%feat) if $debug;

my @col = sort { $header{$b} <=> $header{$a} } keys %header;
print join($sep, @col),"\n" unless $noheader;
for my $id (sort keys %feat) {  # sort is important for Unix "join" command
  print join($sep, map { defined($feat{$id}{$_}) ? $feat{$id}{$_} : $empty } @col),"\n";
}

my $wrote = scalar keys %feat;
print STDERR "Processed $count features\n";
print STDERR "Wrote out $wrote features that had a /$key (set via --key)\n";
print STDERR "Done.\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"format=s",  VAR=>\$format, DEFAULT=>'genbank', DESC=>"Input format"},
    {OPT=>"ftype=s",  VAR=>\$ftype, DEFAULT=>'CDS|RNA', DESC=>"Which feature types"},
    {OPT=>"key=s",  VAR=>\$key, DEFAULT=>'locus_tag', DESC=>"What tag to use primary column"},
    {OPT=>"ignore=s",  VAR=>\$ignore_csv, DEFAULT=>'transl_table,codon_start,translation,db_xref,protein_id,note', DESC=>"Ignore these tags"},
    {OPT=>"keep=s",  VAR=>\$keep_csv, DEFAULT=>'gene,product,EC_number', DESC=>"Include these tags"},
    {OPT=>"noheader!",  VAR=>\$noheader, DEFAULT=>'', DESC=>"Don't print column headings as first line"},
    {OPT=>"empty=s",  VAR=>\$empty, DEFAULT=>'', DESC=>"What to use for empty cells"},
    {OPT=>"boolean=s",  VAR=>\$boolean, DEFAULT=>'yes', DESC=>"What to use for boolean tags when true (eg. /pseudo)"},
    {OPT=>"sep=s",  VAR=>\$sep, DEFAULT=>"\t", DESC=>"Output separator, set to ',' for CSV"},
    {OPT=>"revsuffix=s",  VAR=>\$revsuffix, DEFAULT=>'', DESC=>"Duplicate features and append this to --key eg. 'r' for Nesoni count files"},
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
  print "Usage: $0 [options] file.gbk file2.gbk ... > features.tab\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__ 
locus_tag product inference gene EC_number pseudo anticodon
gene_synonym ribosomal_slippage function ncRNA_class transl_except
