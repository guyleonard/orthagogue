#!/bin/sh 
FILE="/work/mironov/outputs/mpiblast/sprot-2011-11-17/sprot-2011-11-17.blast";
lines=200000;
#lines=500000;
#lines=1000000;
#lines=10000000;
#lines=50000000;
OUT="temp.blast"
head -n $lines $FILE > $OUT;
size=`ls -lh $OUT | awk '{printf $5}'`; # Note: this may change between different architectures.
printf "Found the first %d lines of '%s', and writes them to '%s', givng the size %s\n" $lines $FILE $OUT $size;
head -n 2 $FILE;
./orthaGogue -i temp.blast -s '|' -t 0 -p 1 -c 1 -O result_folder ;