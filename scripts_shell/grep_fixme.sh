#!/bin/sh
STRING="FIXME";
FOLDER=".";
FILE_TYPES_1="*.h";
if [ "$#" -eq 1 ] ; then 
    STRING=$1;
else
    if [ "$#" -eq 2 ] ; then
	STRING=$1;
	FOLDER=$2;
    else
	if [ "$#" -eq 3 ] ; then
	    STRING=$1;
	    FOLDER=$2;
	    FILE_TYPES_1=$3;
	fi
    fi
fi 

echo "Searches in " $FOLDER " of type " $FILE_TYPES_1 " for string: " $STRING
find $FOLDER -name "$FILE_TYPES_1" -exec grep -H $STRING --color -FH {} \; 
if [ "$#" -ne 3 ] ; then
    echo "Searches in " $FOLDER " of type *.cxx for string: " $STRING
    find $FOLDER -name "*.cxx" -exec grep -H $STRING --color -FH {} \; 
fi

exit;
