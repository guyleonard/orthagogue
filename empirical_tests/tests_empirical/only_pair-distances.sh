#!/bin/sh
if [ "$#" -eq 2 ] ; then 
    awk -F " " '{ print $1 "\t" $2 "\t" $12}' $1 1>temp_$1;
    cat temp_$1;
else
    echo "!!  Dot able to run as file not specified."
fi

