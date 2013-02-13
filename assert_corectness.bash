#!/bin/bash
# (1) runs orthAgogue with the correct parameters:
./orthAgogue -i empirical_tests/goodProteins.blast -nn -O empirical_tests/turbo/gp_file/ -b -c 1; 
#./orthAgogue -i empirical_tests/goodProteins.blast -nn -O empirical_tests/turbo/gp_file/ -b -c 1 -pd; 
# (2) Starts the comparison. (If unsure about what happenning, skip the arguments to the script, and the help will provide the answers!): 
perl heuristical_orthology_building/heuristical_orthology_building.pl empirical_tests/goodProteins.blast empirical_tests/turbo/gp_file/ empirical_tests/omcl/gp_file/;

exit;