#!/bin/sh

#
#PBS -A acc-biology
#PBS -lwalltime=4:10:00
#PBS -lnodes=1:ppn=9
#PBS -N speedy_8
module load intelcomp/12.1.0;
module load cmake/2.8.7cmake/2.8.7
cd /home/olekrie/orthaGogue_0.9.9.3/src;

#
# Then the debug version in order to compare time and results:
valgrind --leak-check=full --log-file=valgrind_8.log ./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 8 -O result_folder_debug_8 ;
# wc -l result_folder_debug/*.abc ;
# ls -lh result_folder_debug/*.abc ;
# ./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 2 -O result_folder_debug_2 ;
# wc -l result_folder_debug/*.abc ;
# ls -lh result_folder_debug/*.abc ;
# ./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 3 -O result_folder_debug_3 ;
# wc -l result_folder_debug/*.abc ;
# ls -lh result_folder_debug/*.abc ;
# ./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 4 -O result_folder_debug_4 ;
# wc -l result_folder_debug/*.abc ;
# ls -lh result_folder_debug/*.abc ;
# ./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 5 -O result_folder_debug_5 ;
# wc -l result_folder_debug/*.abc ;
# ls -lh result_folder_debug/*.abc ;