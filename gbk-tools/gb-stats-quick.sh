#!/bin/sh

for FILE in $* ; do

	echo -n "Filename  : "
	echo $FILE

	echo -n "Locus     : "
	head -1 $FILE

	echo -n "Size (bp) : "
	grep ' source ' "$FILE" | sed 's/^.*\.\.//' | head -1
	
	for FEAT in gene CDS tRNA rRNA tmRNA ncRNA misc_RNA ; do
		echo -n "$FEAT      : "
		grep -c "^     $FEAT" $FILE
	done
	
	echo -n "Runs of N : "
	cat "$FILE" | readseq -pipe -i=gb -f=fa 2> /dev/null | \
	  fa-analyze_N_runs.pl - | wc -l
	
	echo
done
