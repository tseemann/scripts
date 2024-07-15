#!/bin/sh

if [ $# -lt 4 -o "$1" = "-h" -o "$1" = "-help" ]; then
  echo "Usage: $0 <reference.gbk> <label> <nesoniDir1 nesoniDir2 ...>"
  exit -1
fi

REF="$1"
shift
echo "Using '$REF' as reference genbank file"

CODE="$1"
shift
echo "Using '$CODE' as file marker"

for A in table compact nexus counts ; do
  echo "### nway $CODE ref $REF as $A on $*"
  nice nesoni nway --require-all no  --gbk "$REF" --as $A $* > nway.pan.${CODE}.$A
  nice nesoni nway --require-all yes --gbk "$REF" --as $A $* > nway.core_diff.${CODE}.$A
  nice nesoni nway --full yes --require-all yes --gbk "$REF" --as $A $* > nway.core_all.${CODE}.$A
done
