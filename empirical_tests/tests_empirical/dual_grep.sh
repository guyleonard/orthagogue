#!/bin/sh

#grep $1 $FILE;
#grep \"$1\" $FILE;
FILE="control_set/averaged_data.abc";
grep $1 $FILE | grep $2;
echo "----------------";
FILE="out_build_o.txt";
grep $1 $FILE | grep $2;
echo "----------------";
FILE="goodProteins.blast"
grep $1 $FILE | grep $2;
echo "----------------.";