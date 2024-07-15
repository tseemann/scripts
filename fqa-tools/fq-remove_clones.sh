#!/bin/bash

# for latest glib
LD_LIBRARY_PATH=/home/sf/devel/C/glib/install/lib 
# real binary
BIN=/home/sf/devel/C/libngs/install/bin/kmers_remove_clonal

if [ "$#" -lt 3 -o "$1" = "-h" -o "$1" = "--help" ]; then
  echo "Usage: $0 <prefix> <left.fq> <right.fq> <options>"
  echo "Options: -c start,end [ -c start2,end2 ... ]"
  exit -1
fi

PREFIX=$1
shift
LEFT=$1
shift
RIGHT=$1
shift
OPTIONS="$*"

if [ ! -r "$LEFT" ]; then
  echo "Unable to read: $LEFT"
  exit -2
fi

if [ ! -r "$RIGHT" ]; then
  echo "Unabled to read: $RIGHT"
  exit -3
fi

if [[ "$LEFT" =~ ".bz2$" ]]; then
  CAT=bzcat
elif [[ "$LEFT" =~ ".gz$" ]]; then
  CAT=zcat
else 
  CAT=cat
fi

echo "Using '$CAT' to read input files."

echo "Launching FIFO $LFIFO : $CAT $LEFT"
LFIFO=seq1.$$
mkfifo $LFIFO
$CAT "$LEFT"  > $LFIFO &

echo "Launching FIFO $RFIFO : $CAT $RIGHT"
RFIFO=seq2.$$
mkfifo $RFIFO
$CAT "$RIGHT" > $RFIFO &

echo "Decloning with options: $OPTIONS"
$BIN -v -o "$PREFIX" "$OPTIONS" $LFIFO $RFIFO
wait

rm -f $LFIFO
rm -f $RFIFO
