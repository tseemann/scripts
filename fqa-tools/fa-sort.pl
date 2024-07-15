#!/usr/bin/env perl
use strict;
use Bio::SeqIO;

my(@Options, $verbose, $by_id, $by_desc, $by_size, $reverse, $numeric, $by_numreads, $by_depth, $by_coverage);
setOptions();

# temp store
my @seq;

# read them in
#print STDERR "Reading...\n";
my $in  = Bio::SeqIO->new(-fh=>\*ARGV, -format=>'Fasta');
while (my $seq = $in->next_seq)  {
  push @seq, $seq;
}

# sort them
if ($by_size) {
  @seq = sort { $a->length <=> $b->length } @seq;
}
elsif ($by_id) {
  @seq = $numeric ? (sort { $a->display_id <=> $b->display_id } @seq)
                  : (sort { $a->display_id cmp $b->display_id } @seq); 
}
elsif ($by_desc) {
  @seq = $numeric ? (sort { $a->desc <=> $b->desc } @seq)
                  : (sort { $a->desc cmp $b->desc } @seq); 
}
elsif ($by_numreads) {
  sub sort_numreads {
    $a->desc =~ m/numreads=(\d+)/; my $aa = $1 || 1;
    $b->desc =~ m/numreads=(\d+)/; my $bb = $1 || 1;
    return $aa <=> $bb;
  }
  @seq = sort { sort_numreads } @seq;
}
elsif ($by_depth) {
  sub sort_depth {
    $a->desc =~ m/numreads=(\d+)/; my $aa = $1 || 1;
    $b->desc =~ m/numreads=(\d+)/; my $bb = $1 || 1;
    return $b->length/$bb <=> $a->length/$aa;
  }
  @seq = sort { sort_depth } @seq;
}
elsif ($by_coverage) {
  # NODE_107_length_51319_cov_108.497887
  sub sort_coverage {
    $a->display_id =~ m/_cov_(\d+\.\d+)/; my $aa = $1 || 1;
    $b->display_id =~ m/_cov_(\d+\.\d+)/; my $bb = $1 || 1;
    return $aa <=> $bb;
  }
  @seq = sort { sort_coverage } @seq;
}
else {
  print STDERR "You must specify one of (--size | --desc | --id | ... etc) to sort by.\n";
  exit 1;
}

# reverse if requested
@seq = reverse @seq if $reverse;

# write output
my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'Fasta');
for my $seq (@seq) {
  $out->write_seq($seq);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"size!",  VAR=>\$by_size, DEFAULT=>0, DESC=>"Sort by size"},
    {OPT=>"id!",  VAR=>\$by_id, DEFAULT=>0, DESC=>"Sort by ID"},
    {OPT=>"desc!",  VAR=>\$by_desc, DEFAULT=>0, DESC=>"Sort by Description"},
    {OPT=>"numreads!",  VAR=>\$by_numreads, DEFAULT=>0, DESC=>"Sort by numreads=XX in desc [454 only]"},
    {OPT=>"depth!",  VAR=>\$by_depth, DEFAULT=>0, DESC=>"Sort by length/numreads in [454 only]"},
    {OPT=>"coverage!",  VAR=>\$by_coverage, DEFAULT=>0, DESC=>"Sort by coverage [Velvet only]"},
    {OPT=>"r|reverse!",  VAR=>\$reverse, DEFAULT=>0, DESC=>"Reverse sort"},
    {OPT=>"n|numeric!",  VAR=>\$numeric, DEFAULT=>0, DESC=>"Numeric sort"},
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
  print "Usage: $0 [options] in.fasta [ more.fasta ... ] > sorted.fasta\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
