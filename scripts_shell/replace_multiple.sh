#!/bin/sh
FIND_STRING="orthaGogue";
REPLACE_STRING="orthAgogue";
ROOT=".";

echo "-\tContinue with FIND_STRING=" $FIND_STRING ", REPLACE_STRING=" $REPLACE_STRING " and ROOT=" $ROOT "? (yes for starting)";
read affirm;
if [ "$affirm" = "yes" ] ; then
    echo "-\tStarts operation:"  $affirm
 else
    echo "--\tAborts the replacements due to argument: " $affirm
    exit;
fi
# Loops through the results from the found-operations
# Remark(1); -F "/" implies that "/" is used as a filed seperator.
# Remark(2): The variable NF in 'awk' is set to the total number of fields in  the input. Usage of this implies that we get the last field, implying the name of the file (and not the path).

for scripts in $(find $ROOT -name "*.h") ; do
  cat $scripts | sed s/$FIND_STRING/$REPLACE_STRING/g > $scripts.tmp;
  cp -f $scripts.tmp $scripts;
  rm -f $scripts.tmp;
done 

for scripts in $(find $ROOT -name "*.cxx") ; do
  cat $scripts | sed s/$FIND_STRING/$REPLACE_STRING/g > $scripts.tmp;
  cp -f $scripts.tmp $scripts;
  rm -f $scripts.tmp;
done

#  echo "cat $scripts | sed s/$FIND_STRING/$REPLACE_STRING/g > $scripts.tmp;" 
#for scripts in $(find . -name "*.h" | awk -F "/" '{print $NF}') ; do