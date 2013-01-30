#------------------------------------------------------------
# REMOVAL OF TEMPORARY FILES:
#------------------------------------------------------------

# Removes the resulting files produced during the oepration:
find . -iname "*.abc" -exec rm -f {} \;
find . -iname "*.mci" -exec rm -f {} \;
find . -iname "*.map" -exec rm -f {} \;
# Removes the debug-files printed for verifying correctness:
find . -iname "out*.txt" -exec rm -f {} \;
find . -iname "err*.txt" -exec rm -f {} \;
find . -iname "temp*.txt" -exec rm -f {} \;
# Removes all blast-files stored in the repository to reduce the size of it:
#mv all.blast all.blast.txt;
#find . -iname "*.blast" -exec rm {} \; 
#mv all.blast.txt all.blast;
# Removes the system-specific cmake-files:
find . -iname "CMakeCache.txt" -exec rm -f {} \;
# Removes the backups of each file:
find . -iname "*~" -exec rm -f {} \;
find . -iname "old*" -exec rm -rf {} \; # If a directory is found, an error is thrown.

