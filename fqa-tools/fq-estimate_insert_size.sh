#!/bin/bash

# load functions common to prokka bash scripts
. `dirname $0`/prokka.functions

if [ $# -ne 3 -o "$1" = "-h" -o "$1" = "-help" ]; then
  usage "refGenome.fna leftReads.fq rightReads.fq"
fi

require_exe maq
require_exe grep
require_exe head

REF=$1
LEFT=$2
RIGHT=$3

if [ ! -r $REF ]; then
  problem "Unable to read refGenome: $REF"
elif [ `head -c 1 $REF` != ">" ]; then
  problem "refSequence must be FASTA format: $REF"
fi

if [ ! -r $LEFT ]; then
  problem "Unable to read leftReads: $LEFT"
elif [ `head -c 1 $LEFT` != "@" ]; then
  problem "leftReads must be Solexa FASTQ format: $LEFT"
fi

if [ ! -r $RIGHT ]; then
  problem "Unable to read rightReads: $RIGHT"
elif [ `head -c 1 $RIGHT` != "@" ]; then
  problem "rightReads must be Solexa FASTQ format: $RIGHT"
fi

require_exe maq

# target directory for everything
DIR=`mktemp -t -d` || problem "Could not make temp directory";
inform "Working dir: $DIR"

# Estimate insert size

SAMPLE=400000
SAMPLELINES=$[ 4 * $SAMPLE ]
inform "Extracting $SAMPLE read pairs ($SAMPLELINES lines) to estimate insert size"
head -$SAMPLELINES $LEFT | maq sol2sanger /dev/stdin $DIR/left.fq
head -$SAMPLELINES $RIGHT | maq sol2sanger /dev/stdin $DIR/right.fq

inform "Converting to MAQ input files"
maq fasta2bfa $REF $DIR/ref.bfa 2> /dev/null
maq fastq2bfq $DIR/left.fq $DIR/left.bfq 2> /dev/null
maq fastq2bfq $DIR/right.fq $DIR/right.bfq 2> /dev/null

inform "Mapping $SAMPLE reads with MAQ"
maq map -a 1000 $DIR/maq.map $DIR/ref.bfa $DIR/left.bfq $DIR/right.bfq 2> $DIR/maq_map.out
#cat $DIR/maq_map.out

INS_MU=`grep 'cal_insert_size' $DIR/maq_map.out | sed 's/^.*: //' | cut -d ' ' -f 1`
INS_MU=`echo $INS_MU | sed 's/\..*$//'`  # round to integer
#INS_MU=$(( $INS_MU - 2*36 ))
INS_SD=`grep 'cal_insert_size' $DIR/maq_map.out | sed 's/^.*: //' | cut -d ' ' -f 3`

# this is to STDERR
inform "Estimated insert size: $INS_MU"
inform "Estimated insert sdev: $INS_SD"

inform "Removing temporary directory/files $DIR"
rm -fr $DIR

inform "Done."

# write to STDOUT too for automatic processing
echo "-ins_length $INS_MU -ins_length_sd $INS_SD"


