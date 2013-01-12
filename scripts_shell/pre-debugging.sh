#!/bin/sh

echo""
echo "Runs the script enabling simple debugging, setting the parameters in the CMakeLists.txt-files" ; 

# Removes the system-specific cmake-files:
find . -iname "CMakeCache.txt" -exec rm {} \;

#------------------------------------------------------------
# STATUS MESSAGES (Hopefully none file should be found.)
#------------------------------------------------------------

# Finds the files above 10MB, producing a warning. 
find . -type f -size +50000k -exec ls -lh {} \; | awk '{ print "\tConsider removing " $7 ": " $5 }' 
find -name "*.cxx" -o -name "*.h" -exec grep -H "FIXME" --color -FH {} \; 
#find -name "*.cxx" -o -name "*.h" -exec grep -H "TODO" --color -FH {} \; 
find . -name "*.cxx" -exec wc -l {} \; | awk '{ SUM += $1 } END { print "Line-count for the *.cxx-files is " SUM " lines."}'
find . -name "*.h" -exec wc -l {} \; | awk '{ SUM += $1 } END { print "Line-count for the *.h-files is " SUM " lines."}'

#cmake .; make;