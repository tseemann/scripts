package FastQ;

use Exporter qw(import);
@EXPORT_OK = qw(read_fastq write_fastq write_fasta read_fasta
                write_qual int_qual int_to_qual reverse_complement
		trim_fastq assert trim_fastq_NEW);

use Data::Dumper;
use Carp;

#----------------------------------------------------------------------

sub assert {
  my($condition, $errormsg) = @_;
  unless ($condition) {
    print STDERR "ERROR: $errormsg\n";
    exit 1;
  }
}

#-------------------------------------------------------------------------

sub read_fastq {
  my($fh) = @_;
#  return if eof $fh;
  my @s = map { scalar(<$fh>) } 1..4;   # read 4 lines
  chomp(@s);                            # remove \n from ends
  $s[0] = substr($s[0], 1);             # remove @ from ID
  splice @s,2,1;                        # remove 3rd line
  return \@s;                           # return 3-tuple [ ID,seq,qual ]
}

#-------------------------------------------------------------------------

sub read_fasta {
  my($fh) = @_;
#  return if eof $fh;
  my @s = map { scalar(<$fh>) } 1..2;   # read 2 lines
  chomp(@s);                            # remove \n from ends
  $s[0] = substr($s[0], 1);             # remove > from ID
  $s[2] = 'B'x(length($s[1]));          # fake a 3rd quality line?
  return \@s;                           # return 3-tuple [ ID,seq,qual ]
}

#-------------------------------------------------------------------------

sub write_fastq {
  my($fh,$s) = @_;
  print $fh '@',$s->[0],"\n",$s->[1],"\n+\n",$s->[2],"\n";
}

#-------------------------------------------------------------------------

sub write_fasta {
  my($fh,$s) = @_;
  print $fh '>',$s->[0],"\n",$s->[1],"\n";
}

#-------------------------------------------------------------------------

sub int_qual {
  my($offset, $sobj) = @_;
  return map { ord($_)-$offset } (split m//, $sobj->[2]);
}

#-------------------------------------------------------------------------

sub int_to_qual {
  my($offset, $array) = @_;
  return join '', map { chr($_+$offset) } @$array;
}

#-------------------------------------------------------------------------

sub write_qual {
  my($fh,$s) = @_;
  my @q = int_qual($s);
  print $fh '>',$s->[0],"\n@q\n";
}

#  for (int i=0; i < seq->len; i++) {
#    seq->seq[i]  = COMPLEMENT( old.seq[seq->len - i - 1] );
#    seq->qual[i] = old.qual[seq->len - i - 1];  // quality does not 'complement'
#  }

#-------------------------------------------------------------------------

sub reverse_complement {
  my($s) = @_;
#  print Dumper($s);
  # 0=id, 1=seq, 2=qual
  $s->[1] = reverse $s->[1];
  $s->[1] =~ tr/ATGCatgc/TACGtacg/;
  $s->[2] = reverse $s->[2];
#  print Dumper($s); exit;
  return $s;
}

#-------------------------------------------------------------------------

sub trim_fastq_NEW {
  my($s, $minq) = @_;

  croak "need minq parameter" unless $minq;
  my $L = length($s->[1]);

  my @seq;
  my $begin = 0;
  
  while ($begin < $L) 
  {
    while ($begin < $L and ord(substr($s->[2],$begin,1))-64 < $minq) {
      $begin++;
    }
    last if $begin >= $L;

    my $end = $begin;
    while ($end < $L and ord(substr($s->[2],$end,1))-64 >= $minq) {
      $end++;
    }
    $end--;  # we always go one over?

    push @seq, [ 
      $s->[0], 
      substr($s->[1],$begin,$end-$begin+1),
      substr($s->[2],$begin,$end-$begin+1),
    ];
    
    $begin = $end+1;
  }
  return @seq;
}  

#-------------------------------------------------------------------------

sub trim_fastq {
  my($s, $minq) = @_;

  croak "need minq parameter" unless $minq;
  my $L = length($s->[1]);

  my $begin = 0;
  while ($begin < $L and ord(substr($s->[2],$begin,1))-64 < $minq) {
    $begin++;
  }
  return if $begin >= $L;
  
  my $end = $begin;
  while ($end < $L and ord(substr($s->[2],$end,1))-64 >= $minq) {
    $end++;
  }
  $end--;  # we always go one over?
  
  return [ 
    $s->[0], 
    substr($s->[1],$begin,$end-$begin+1),
    substr($s->[2],$begin,$end-$begin+1),
  ];
}  

#-------------------------------------------------------------------------

1;
