#!/bin/sh
echo ""
make; # As a developer I always have to compile in order to ensure the latest changes are in effect.
echo "Compares the output: the line numbers should equal."
./orthaGogue -i all_vm.blast -c 1 1>enkel.txt;  # The default running
mv all.abc all_default.abc; # Stores the result for later comparison for the default run.
./orthaGogue -i all_vm.blast -c 1 -nii 1>nii_enkelt.txt; # Runs with restricted defenition
mv all.abc all_nii.abc; # Stores the result for later comparison for the default run.
diff nii_enkelt.txt enkel.txt; # Compares them: no output is given if they are equal.
wc -l all_default.abc all_nii.abc; # Looks upon the line numbers: the "all_default.abc" should by def be greater- or equal the "all_nii.abc" file.