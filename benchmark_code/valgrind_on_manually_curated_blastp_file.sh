#!/bin/sh
./tests_empirical/produce_smaller_file_for_verific.sh;
valgrind --leak-check=full --log-file=valgrind.log ./orthaGogue -i temp.blast -s '|' -t 0 -p 1 -c 1 -O result_folder;