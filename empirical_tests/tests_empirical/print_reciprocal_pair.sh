#!/bin/sh
# Script written for verifying the reciprocality.
if [ "$#" -eq 4 ] ; then 
    # Note: In order to remove ducplicates and containing alle the values, 'sort' is used with parameters for the columns of id with the additional requirement of uniqueness.
    #cmd="grep \"$1\" $3 | grep \"$2\" | sort --key=1,2 -u | awk -F \" \" '{ print $1 \"\t\" $2 \"\t\" $12}'"; echo "$cmd"
    grep "$1" $3 | grep "$2" | sort --key=1,2 -u | awk -F " " '{ print "(blastp)  " $1 "\t" $2 "\t" $12}';# #> $fname;
    grep "$1" $4 | grep "$2"|  awk -F " " '{ print "(turbOO)  " $1 "\t" $2 "\t" $3}';# #> $fname;
    echo "-------------------------------------------------------"
else
    echo "!!  Not able to run as file not specified."
fi
