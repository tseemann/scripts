#!/bin/sh

if [ $# -lt 1 ]; then
  echo "Usage: $0 <velvet_dir>" 
  exit
fi

DIR=$1

if [ ! -d "$DIR" ]; then
  echo "$DIR is not a directory"
  exit
fi

if [ ! -r "$DIR/Log" ]; then
  echo "$DIR does not contain a Log file"
  exit
fi

echo
echo "VERSION"
echo
velveth | egrep '^(Version|CATEGOR|MAXKMER)'

echo
echo "LOG"
echo
cat $DIR/Log

echo
echo "BACKTRACE"
echo
echo "[velveth]"
gdb velveth core.* << EOF | grep '^#'
bt
EOF
echo "[velvetg]"
gdb velvetg core.* << EOF | grep '^#'
bt
EOF

echo
echo "EXECUTABLES"
echo
file `which velvetg`
file `which velveth`

echo
echo "SYSTEM"
echo
uname -a
free -o | grep ':'
CORES=`grep -c '^processor' /proc/cpuinfo`
echo -n "$CORES x "
grep '^model name' /proc/cpuinfo | head -1


echo



