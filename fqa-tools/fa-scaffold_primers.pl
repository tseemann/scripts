#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Tools::Run::Primer3;
use Data::Dumper;
use File::Temp;
use Fatal;
use Bio::Seq::PrimedSeq;

my(@Options, $verbose, $prefix, $dist, $length, $zone, $mispriming);
setOptions();

my $in = Bio::SeqIO->new(-file=>$ARGV[0], -format=>'Fasta');
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');

my $dna = '';
my $pos = 0;
my @primer;
my $right;
my $left;

while (my $seq = $in->next_seq) {
  my $id = $seq->display_id;
  print STDERR "Loaded: $id\n";
  my @s = split m/N+/, uc($seq->seq);
  print STDERR "Regions: ", scalar @s, "\n";
  my $count=0;
  for my $s (@s) {
    $dna .= $s;
    my $L = length($s);
    print STDERR "Contig $count: $L bp\n" if $verbose;
    
    if ($count > 0) {
      $right = design_primer("$id.$count.R", "RIGHT",  $pos+$dist) unless $count==0;
      if ($right and $left) {
        $out->write_seq($left, $right);
      }
      else {
        print STDERR "WARNING: could not design pair near contig $count pos $pos\n";
      }	
    }
    $count++;
    $left = design_primer("$id.$count.L", "LEFT", $pos+$L-$dist-$length-$zone);

    $pos += $L;
  }
} 

#---------------------------------------------------------------------

sub design_primer {
  my($name, $dir, $pos) = @_;
  my $p3 = Bio::Tools::Run::Primer3->new;
  $p3->executable('/bio/sw/primer3-1.1.4/primer3_core');
#  my $args = $p3->arguments;
#  foreach my $key (keys %{$args}) {print "$key\t", $args->{$key}, "\n"}
#  exit;
  $p3->add_targets(
    'PRIMER_SEQUENCE_ID'  => $name,
    'SEQUENCE'            => $dna,
    'INCLUDED_REGION'     => "$pos,$zone",
    'PRIMER_TASK'         => lc("pick_${dir}_only"),
    'PRIMER_OPT_SIZE'     => $length,
    'PRIMER_EXPLAIN_FLAG' => 1,
  );
  $p3->add_targets('PRIMER_MAX_TEMPLATE_MISPRIMING', $mispriming) if $mispriming;
  print STDERR "Running: $name $dir $pos ...\n" if $verbose;
  my $res = $p3->run;
#  print Dumper($res);
  print STDERR "Primers found: ", $res->number_of_results, "\n" if $verbose;
  for my $index (0 .. $res->number_of_results-1) {
#    print STDERR "Testing: #$index\n";
    my $p = $res->primer_results($index);
#    print Dumper($p); exit;
    my $ps = Bio::Seq->new(-id=>$name,-seq=>$p->{"PRIMER_".uc($dir)."_SEQUENCE"});
    return $ps if is_unique($ps, $dir);
  }
}

#---------------------------------------------------------------------

sub is_unique {
  my($primer, $dir) = @_;
  print STDERR "Checking ", $primer->display_id, " ",$primer->seq, " against $ARGV[0]\n";
  my $fh = File::Temp->new(SUFFIX=>'.fasta');
  print $fh ">".$primer->display_id."\n".$primer->seq."\n";
  my $exo = Bio::SearchIO->new(
    -format=>'exonerate',
    -file=>"exonerate --model affine:local $fh $ARGV[0] |",
  );
  my $res = $exo->next_result;
#  warn "bad exo result for ".$primer->display_id if not defined $res;
  return unless $res;
#  print Dumper($res); exit;
  printf STDERR "Found %d hits to template\n", $res->num_hits if $verbose;

  die if $res->num_hits > 1;

  if ($res->num_hits == 1) {
    my $hit = $res->next_hit;
#    die "bad strand for left!" if $dir =~ m/left/i and $hit->strand('hit') < 0;
#    die "bad strand for right!" if $dir =~ m/right/i and $hit->strand('hit') > 0;
    return if $dir =~ m/left/i and $hit->strand('hit') < 0;
    return if $dir =~ m/right/i and $hit->strand('hit') > 0;
    my $hsp = $hit->next_hsp;
    $primer->desc( join ' ', 
      'unique hit',
      $primer->length,
      'bp',
      'on',
#      (100*$hsp->frac_identical('hsp')).'%id',
      $hit->name, 
      $hit->start('hit'), 
      $hit->end('hit'), 
      $hit->strand('hit'),
    );
    return $primer;
  }
  print STDERR "Bummer, ", $res->num_hits, " hits found!\n";
  return;
}

#---------------------------------------------------------------------

sub runcmd {
  my @cmd = @_;
  print STDERR "Running: @cmd\n";
  system(@cmd)==0 or die "ERROR RUNNING: @cmd";
}  

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose"},
    {OPT=>"p|prefix=s",  VAR=>\$prefix, DEFAULT=>'out', DESC=>"Output prefix"},
    {OPT=>"d|dist=i",  VAR=>\$dist, DEFAULT=>100, DESC=>"Distance from N regions"},
    {OPT=>"l|length=i",  VAR=>\$length, DEFAULT=>20, DESC=>"Optimal primer length"},
    {OPT=>"z|zone=i",  VAR=>\$zone, DEFAULT=>100, DESC=>"Primer allowed zone"},
    {OPT=>"m|mispriming=f",  VAR=>\$mispriming, DEFAULT=>0, DESC=>"PRIMER_MAX_TEMPLATE_MISPRIMING (zero = NONE)"},
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
  print "Usage: $0 [options] [-p prefix] scaffolds.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------

__DATA__

SEQUENCE_ID=example
SEQUENCE_TEMPLATE=GTAGTCAGTAGACNATGACNACTGACGATGCAGACNACACACACACACACAGCACACAGGTATTAGTGGGCCATTCGATCCCGACCCAAATCGATAGCTACGATGACG
SEQUENCE_TARGET=37,21
PRIMER_TASK=pick_detection_primers
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=75-100
P3_FILE_FLAG=1
SEQUENCE_INTERNAL_EXCLUDED_REGION=37,21
PRIMER_EXPLAIN_FLAG=1
=
