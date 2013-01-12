#!/bin/bash
make 1>out_make.txt; mpirun -np 2 ./orthAgogue -i empirical_tests/goodProteins.blast -A -dbs 10024 -c 2
exit;