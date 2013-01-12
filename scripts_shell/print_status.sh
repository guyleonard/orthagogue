#------------------------------------------------------------
# STATUS MESSAGES (Hopefully none file should be found.)
#------------------------------------------------------------

# Finds the files above 10MB, producing a warning. 
find . -type f -size +5000k -exec ls -lh {} \; | awk '{ print "\tConsider removing " $9 ": " $5 }' 
find . -type f -size +5000k -exec ls -lh {} \; | awk '{ printf $9 " " }' 
echo ""
#find . -type f -size +50000k -exec ls -lh {} \; | awk '{ print "\tConsider removing " $7 ": " $5 }' 
# Locates the FIXME's in the files:
find -name "*.cxx" -exec grep -H "FIXME" --color -FH {} \; 
find -name "*.h" -exec grep -H "FIXME" --color -FH {} \; 
find -name "*.txt" -exec grep -H "FIXME" --color -FH {} \; 
find -name "*.cmake" -exec grep -H "FIXME" --color -FH {} \; 

# Sums the lines written:
find . -name "*.cxx" -exec wc -l {} \; | awk '{ SUM += $1 } END { print "Line-count for the *.cxx-files is " SUM " lines."}'
find . -name "*.h" -exec wc -l {} \; | awk '{ SUM += $1 } END { print "Line-count for the *.h-files is " SUM " lines."}'

