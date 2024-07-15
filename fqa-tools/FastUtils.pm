package FastUtils;

use Exporter qw(import);
@EXPORT_OK = qw(read_one_seq read_fastq read_fastq_Bio read_fasta read_fasta_Bio
				write_fastq write_fasta write_qual int_qual int_to_qual
				trim_fastq assert trim_fastq_NEW detect_filetype);

use Data::Dumper;
use Carp;
use Bio::Seq;

#----------------------------------------------------------------------

sub assert {
  my($condition, $errormsg) = @_;
  unless ($condition) {
    print STDERR "ERROR: $errormsg\n";
    exit 1;
  }
}

#-------------------------------------------------------------------------
# simple sub to try and detect a seqeunce file type from:
# (single line fasta, single line fastq, multiline fasta, multiline fastq)
sub detect_filetype {
	my($f) = @_;
	open $fh, "<", $f;
	my @s = map{ scalar (<$fh>) } 1..4;
	close $fh;
	my $filetype = "";
	if($s[0] =~ /^\@/){
		$filetype = "Fastq";
	}
	elsif($s[0] =~ /^>/){
		$filetype = "Fasta";
	}
	else {
		$filetype = "Unknown";
	}
	unless($filetype eq "Unknown" || $s[2] =~ /^[>\+]/){
		$filetype .= "_Bio";
	}
	return $filetype;
}

#-------------------------------------------------------------------------
# picks the correct read_fast? sub depending on $bio and $fasta...
# but needs to do this for every seq...  TODO: Dynamically set up correct
# read subroutine once when deciding on .................
sub read_one_seq {
	my ($fh, $bio, $fasta) = @_;
	if($bio){
		$fasta ? return read_fasta_Bio($fh) : return read_fastq_Bio($fh);
	}
	else {
		$fasta ? return read_fasta($fh) : return read_fastq($fh);
	}
	return 0;
}

#-------------------------------------------------------------------------
# reads a fastq entry in a fastq file assuming only one line per sequence
sub read_fastq {
  my($fh) = @_;
  my @s = map { scalar(<$fh>) } 1..4;   # read 4 lines
  chomp(@s);                            # remove \n from ends
  $s[0] = substr($s[0], 1);             # remove @ from ID
  splice @s,2,1;                        # remove 3rd line
  return \@s;                           # return 3-tuple [ ID,seq,qual ]
}

#-------------------------------------------------------------------------
# reads a fastq entry in a fastq file assuming its given a Bio::Seq object
sub read_fastq_Bio {
  my($sio) = @_;
  my $seq = <$sio>;
  return unless defined $seq;
  my @s = ($seq->id(), $seq->seq(), $seq->qual()->qual_text());
  return \@s;
}

#-------------------------------------------------------------------------
# reads a fasta entry in a fasta file assuming only one line per sequence
sub read_fasta {
	my($fh) = @_;
	my @s = map { scalar(<$fh>) } 1..2;   # read 2 lines
	chomp(@s);                            # remove \n from ends
	$s[0] = substr($s[0], 1);             # remove > from ID
	return \@s;                           # return 2-tuple [ ID,seq ]
}

#-------------------------------------------------------------------------
# reads a fasta entry in a fasta file assuming its given a Bio::Seq object
sub read_fasta_Bio {
  my($sio) = @_;
  my $seq = <$sio>;
  return unless defined $seq;
  my @s = ($seq->id(), $seq->seq());
  return \@s;
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
  my($s) = @_;
  return map { ord($_)-64 } (split m//, $s->[2]);
}

#-------------------------------------------------------------------------

sub int_to_qual {
  return join '', map { chr($_+64) } @_;
}

#-------------------------------------------------------------------------

sub write_qual {
  my($fh,$s) = @_;
  my @q = int_qual($s);
  print $fh '>',$s->[0],"\n@q\n";
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
