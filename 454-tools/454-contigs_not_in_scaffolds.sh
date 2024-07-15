#!/bin/sh

if [ ! -r "$1/454Scaffolds.txt" ]; then
	echo "Synopsis: Writes '454ContigsNotInScaffolds.fna' FASTA - those contigs not present in scaffolds"
	echo "Usage: $0 <454_assembly_dir>"
	exit
fi

if [ ! -r "$1/454AllContigs.fna" ]; then
	echo "Could not see '454AllContigs.fna' in folder '$1'"
	exit
fi


DEST="$1/454ContigsNotInScaffolds.fna"
grep contig "$1/454Scaffolds.txt" | cut -f 6 | fa-extract.pl -i "$1/454AllContigs.fna" -l /dev/stdin --false $DEST 2> /dev/null
echo "Wrote: $DEST"
echo -n "Contigs: "
grep -c '>' $DEST

SALL="$1/454All.fna"
cat $1/454Scaffolds.fna $1/454ContigsNotInScaffolds.fna | fa-filter-size.pl --min 0 - 1> $SALL 2> /dev/null
echo "Wrote: $SALL"
echo -n "Contigs: "
grep -c '>' $SALL

SALL="$1/454All250.fna"
cat $1/454Scaffolds.fna $1/454ContigsNotInScaffolds.fna | fa-filter-size.pl --min 250 - 1> $SALL 2> /dev/null
echo "Wrote: $SALL"
echo -n "Contigs: "
grep -c '>' $SALL
