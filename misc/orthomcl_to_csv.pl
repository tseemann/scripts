#!/usr/bin/env perl
use strict;
use List::Util qw(max first);
use List::MoreUtils qw(uniq);
use Data::Dumper;
#use Math::Combinatorics;
#use Algorithm::Combinatorics qw(:all);
use List::PowerSet qw(powerset);
use Text::CSV;
use Bio::SeqIO;

my(@Options, $debug);
setOptions();

my $dir = shift @ARGV;
-d $dir or usage();

#ORTHOMCL7(7 genes,3 taxa):       Sal1296_cds05921(1) Sal1296_cds05922(1) Sal1309_cds05726(2) Sal1309_cds05727(2) Sal1309_cds05728(2) Sal1635_cds05661(3) Sal1635_cds05662(3)
#ORTHOMCL8(6 genes,3 taxa):       Sal1296_cds00053(1) Sal1296_cds00786(1) Sal1309_cds00953(2) Sal1309_cds06374(2) Sal1635_cds00051(3) Sal1635_cds00736(3)

my $ABSENT = 'absent';
my $T = 0;
my $N = 0;

my %group_of;
my %ids_of;
my %patt_of;
my %genome;

sub numeric { $a <=> $b };

open ORTHO, "$dir/all_orthomcl.out";

while (<ORTHO>) {
  chomp;
  next unless m/^ORTHOMCL(\d+)\((\d+)\s+genes,(\d+)\s+taxa\):\s+(.*)$/i;
  my($code,$genes,$taxa,$list) = ($1,$2,$3,$4);
  my @id = split m/\s+/, $list;
  die "mismatch on gene counts: $genes <> @id" if @id != $genes;
  $T = max($taxa, $T); # deduce no. of genomes
  $N++; # no. of groups
  print "$code @id\n" if $debug;
  map { $group_of{$_} = $code } @id;
  $ids_of{$code} = [ @id ];
  my @patt;
  foreach (@id) {
    m/\((\w+)\)$/;
    push @patt, $1;
    $genome{$1}++;
  }
  $patt_of{$code} = set_pattern(@patt);
  print "$code | ",$patt_of{$code}," | @patt\n" if $debug;
}

die "mismatch on group counts: $N <> ?" if $N != scalar(keys %ids_of);
#print Dumper(\%group_of, \%ids_of, \%patt_of, \%genome) if $debug;
#print Dumper(\%genome); 
#print "N=$N\n";
#print Dumper(\%ids_of); 
#print Dumper(\%group_of); 
#print Dumper(\%patt_of); 
#exit;


# GO BACK AND INCLUDE EXCLUDED SINGLETONS ?

open GGFILE, "$dir/tmp/all.gg";
while (<GGFILE>) {
  chomp;
  my($gid,@id) = split ' ';
  $gid =~ s/:$//;
  for my $id (@id) {
    my $oid = "$id($gid)";
    next if exists $group_of{$oid};
    $group_of{$oid} = ++$N;
    $patt_of{$N} = $gid;
    $ids_of{$N} = [ $oid ];
    $genome{$gid}++;
  }
}
#print Dumper(\%genome); 
#print "N=$N\n";
#print Dumper(\%ids_of); 
#print Dumper(\%patt_of); 


my %product_of;
my %length_of;
my %num_x_in;
my $fasta = "$dir/tmp/all.fa";
print STDERR "Opening FASTA: $fasta\n";
my $fin = Bio::SeqIO->new(-file=>$fasta, -format=>'Fasta');
while (my $seq = $fin->next_seq) {
  my $id = $seq->display_id;
  $product_of{ $id } ||= $seq->desc;
  $length_of{ $id } ||= $seq->length;
  $num_x_in{ $id } ||= ($seq->seq =~ tr/Xx/Xx/);
}
#print STDERR Dumper(\%product_of);

my $csv = Text::CSV->new;

my @gid = sort keys %genome;
my $x = powerset(1 .. $T);
my @set = reverse sort { scalar(@$b) <=> scalar(@$a) } @$x;
#print Dumper(\%patt_of);

sub gene_sort {
  $ids_of{$a}->[0] cmp $ids_of{$b}->[0]
}

for my $set (@set) {
  next unless @$set > 0;
  print STDERR Dumper($set) if $debug;
  my @setstr = map { $gid[$_-1] } @$set;
  my $set_patt = set_pattern(@setstr);
  print "\n";
#  print_csv("SUBSET", "GROUP", @gid, 'LENGTH', 'NUM_X', 'PRODUCT');
  print_csv("COUNTER", "SUBSET", @gid, 'LENGTH', 'NUM_X', 'PRODUCT');
  my $counter=0;
  for my $g (sort gene_sort keys %patt_of) {
    next unless $set_patt eq $patt_of{$g};
    my @ids = partition_ids($T, @{$ids_of{$g}});
    my($id) = split m/\s+/, (first { $_ ne $ABSENT } @ids);
    my $prod = $product_of{$id} || 'unknown';
    my $len = $length_of{$id} || 'unknown';
    my $numx = $num_x_in{$id} || '0';
    $prod =~ s/\s+\(.*$//;
#    print_csv($set_patt, $g, @ids, $prod);
    print_csv(++$counter, $set_patt, @ids, $len, $numx, $prod);
  }
}

sub print_csv {
  return unless @_;
  $csv->combine(@_) or die $csv->error_diag;
  print $csv->string,"\n";
}

sub partition_ids {
  my($T,@id) = @_;
  my @res;
  for my $g (sort keys %genome) {
    my $r = (join " ", grep { m/\($g\)$/ } @id) || $ABSENT;
    $r =~ s/\(\w+\)//g;
    push @res, $r;
  }
  return @res;
}

sub set_pattern {
  join('+', uniq sort @_);
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
#    {OPT=>"fasta=s",  VAR=>\$fasta, DEFAULT=>'', DESC=>"Fasta file"},
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
  print "Usage: $0 [options] orthomcl_result_directory > orthologs.csv\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
