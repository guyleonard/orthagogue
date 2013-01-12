#!/bin/sh

# Finds the files above 10MB:
find . -type f -size +50000k -exec ls -lh {} \; | awk '{ print $9 ": " $5 }' 

# Greps the things to do before release, with highest priority to the FIXME's
grep -H -r "FIXME"  *.h
grep -H -r "FIXME"  *.cxx
grep -H -r "FIXME"  *.txt *.cmake *.sh *.bash
#grep -H -r "TODO"  .