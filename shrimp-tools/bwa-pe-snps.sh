#!/bin/sh

EXENAME=`basename $0`
NUMPARAM=4
CPUS=8

if [ $# -lt $NUMPARAM -o "$1" = "-h" -o "$1" = "--help" ]; then
  cat << HELPTEXT 1>&2
NAME
  $EXENAME - Find SNPs using BWA and SAMtools
SYNOPSIS
  $EXENAME reference.fna left.fq right.fq folder
PARAMETERS
  reference.fna  FASTA file containing one or more reference sequences
  left.fq        FASTQ file for left mates
  right.fq       FASTQ file for right mates
  folder         Destination folder for all results and intermediate files
OUTPUT
  Look in folder/ for the goodies
AUTHOR
  Torsten Seemann <torsten.seemann@infotech.monash.edu.au>
HELPTEXT
  exit -1
fi

REF=$1
LEFT=$2
RIGHT=$3
DEST=$4

if [ -d "$DEST" ]; then
  echo "ERROR: Destination folder $DEST already exists."
  exit -2
fi

echo "Creating folder: $DEST"
mkdir -p "$DEST"

# bwa index

cp -f "$REF" "$DEST/ref.fna"
pushd "$DEST"
bwa index ref.fna
popd

ls -lsa "$DEST"

# bwa aln

PREFIX="$DEST/ref.fna"
bwa aln -t $CPUS -f "$DEST/left.sai" "$PREFIX" "$LEFT"
bwa aln -t $CPUS -f "$DEST/right.sai" "$PREFIX" "$RIGHT"

# bwa sampe

bwa sampe -P -f "$PREFIX.sam" "$PREFIX" \
	"$DEST/left.sai" "$DEST/right.sai" \
	"$LEFT" "$RIGHT"

# convert sam to bam

samtools view -b -T "$PREFIX" -o "$PREFIX.bam.unsorted" "$PREFIX.sam"

# sort the bam

samtools sort "$PREFIX.bam.unsorted" "$PREFIX"

# raw SNPs

samtools pileup -c -v -N 1 -f "$PREFIX" -N 1 "$PREFIX.bam" > "$PREFIX.pileup"

# filtered SNPs

samtools.pl varFilter "$PREFIX.pileup" | awk '$8>=10' > "$PREFIX.snps"
grep '\*' "$PREFIX.snps" > "$PREFIX.indels"

# Cleanup

ls -lsa "$PREFIX.snps" "$PREFIX.indels"

