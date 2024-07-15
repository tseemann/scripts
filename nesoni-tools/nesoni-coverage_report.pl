#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::SeqIO;
use File::Spec;

my(@Options, $verbose, $mincov);
setOptions();
$|=1; # autoflush stdout

for my $dir (@ARGV) {
  -d $dir and -r "$dir/report.txt" or die "'$dir' is not a valid Nesoni folder";
  my $folder = File::Spec->rel2abs($dir);
  printf STDERR "Nesoni folder: $folder\n", 
  #gi_57650036_ref_NC_002951_2_-ambiguous-depth.userplot
  my $suffix = "-ambiguous-depth.userplot";
  for my $plot (glob("$dir/*$suffix")) {
    print STDERR "Processing: $plot\n";
    print "FOLDER\n\t$folder\n";
    my $ref = substr($plot, 0, -length($suffix));
    printf "REFERENCE\n\t%s\n", $ref;
    my $bp = lines($plot);
    printf "\t$bp bp\n";
  #  print "READS\n\t"; 
  #  cat("$dir/config.txt");
    print "DEPTH\n";
    for my $func (qw(mean median min max)) {
      print "\t$func: "; 
      print qx(math $func $plot);
    }
    print "COVERAGE\n";
    for my $m (sort { $a<=>$b} (1, 5, 10, 20, 50, 100, 200)) {
      my $c = qx(math above $m $plot | math count);
      chomp $c;
      printf "\tcovered by depth >= $m : %d bp (%.2f%%)\n", $c, ($c*100/$bp);
    }
#    print "ALIGNMENT\n";
#    cat("$dir/shrimp_log.txt");
  }
}

sub cat {
  
  for (@_) {
    system("cat", $_);
  }
}

sub lines {
  my($n) = qx(cat \Q$_[0]\E | wc -l);
  chomp $n;
  return $n;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"c|mincov=f",  VAR=>\$mincov, DEFAULT=>0, DESC=>"Minimum coverage to consider 'covered'"},
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
  print "Usage: $0 [options] <nesoni_dir>\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
