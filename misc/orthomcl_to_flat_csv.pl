#!/usr/bin/env perl
use strict;
use List::Util qw(max sum first);
use List::MoreUtils qw(uniq);
use Data::Dumper;
#use Math::Combinatorics;
#use Algorithm::Combinatorics qw(:all);
#use List::PowerSet qw(powerset);
use Text::CSV;
use Bio::SeqIO;

my(@Options, $debug);
setOptions();
$|=1;

my $dir = shift @ARGV;
-d $dir or usage();

#ORTHOMCL7(7 genes,3 taxa):       Sal1296_cds05921(1) Sal1296_cds05922(1) Sal1309_cds05726(2) Sal1309_cds05727(2) Sal1309_cds05728(2) Sal1635_cds05661(3) Sal1635_cds05662(3)
#ORTHOMCL8(6 genes,3 taxa):       Sal1296_cds00053(1) Sal1296_cds00786(1) Sal1309_cds00953(2) Sal1309_cds06374(2) Sal1635_cds00051(3) Sal1635_cds00736(3)

my $ABSENT = 'absent';
my $T = 0;
my $N = 0;

my %group_of;
my %ids_of;
my %gid_of;
my %gid_count;

#
# READ IN ORTHOLOG ASSIGNMENTS
#

open ORTHO, "$dir/all_orthomcl.out";
while (<ORTHO>) {
  chomp;
  next unless m/^ORTHOMCL(\d+)\((\d+)\s+genes,(\d+)\s+taxa\):\s+(.*)$/i;
  my($code,$genes,$taxa,$list) = ($1,$2,$3,$4);
  my @id;
  for my $s (split m/\s+/, $list) {
    $s =~ m/^(.*?)\((.*)\)$/ or die "invalid sequence name: $s";
    push @id, $1;
    $gid_of{$1} = $2;
    $gid_count{$2}++;
    $group_of{$1} = $code;
    push @{$ids_of{$code}}, $1;
  }
  die "mismatch on gene counts: $genes <> @id" if @id != $genes;
  print "$code @id\n" if $debug;
}
$N = scalar keys %ids_of;

#print Dumper(\%gid_count);print Dumper(\%gid_of);

#
# GO BACK AND INCLUDE THE EXCLUDED SINGLETONS
#

open GGFILE, "$dir/tmp/all.gg";
while (<GGFILE>) {
  my $line=$_;
  chomp;
  my($gid,@id) = split ' ';
  $gid =~ s/:$//;
  for my $id (@id) {
#    print STDERR "checking: $id\n";
    # not singleton so taken care of already
    next if exists $group_of{$id}; 
    # make new group with one member
    $N++;
#    print STDERR "singleton: $N $gid $id\n";
    $group_of{$id} = $N;
    $gid_of{$id} = $gid;
    $ids_of{$N} = [ $id ];
    $gid_count{$gid}++;
  }
}
$T = scalar keys %gid_count;

printf STDERR "Read %d genes from %d genomes in %d ortholog groups.\n", 
	scalar(keys %group_of), $T, $N;

#print Dumper(\%gid_count, \%group_of, \%ids_of, \%gid_of); exit;

#
# READ IN FASTA FILE TO GET GENE PRODUCTS AND LENGTHS
#

my %seq;
my $fasta = "$dir/tmp/all.fa";
my $fin = Bio::SeqIO->new(-file=>$fasta, -format=>'Fasta');
while (my $seq = $fin->next_seq) {
  print STDERR "\rReading gene sequences: $." if $. % 500 == 0;
  $seq{ $seq->display_id } = $seq;
}
print STDERR "\n";

#
# OUTPUT A NICE FLAT TABLE
#

my $csv = Text::CSV->new;
my @gid = sort keys %gid_count;

my $counter=0;
print_csv('ortholog', 'product', 'length', @gid);
for my $grp (sort { $a <=> $b } keys %ids_of) {
  my $prod = '';
  my $len = 0;
  my @has;
  for my $gid (@gid) {
    my @members = @{ $ids_of{$grp} };
    @members > 0 or die "group $gid has no members!";
    $prod ||= $seq{$members[0]}->desc;
    $len ||= $seq{$members[0]}->length;
    my @rep = grep { $gid_of{$_} eq $gid } @members;
    my $rep = @rep > 0 ? join(',',@rep) : '-';
    push @has, $rep; 
  }
  print_csv(++$grp, $prod, $len, @has);
}
print STDERR "Done.\n";

sub print_csv {
  return unless @_;
  $csv->combine(@_) or die $csv->error_diag;
  print $csv->string,"\n";
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
