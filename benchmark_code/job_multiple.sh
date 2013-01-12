#!/bin/sh

#
#PBS -m abe
#PBS -A acc-biology
#PBS -lwalltime=0:40:00
#PBS -lnodes=1:ppn=7
#PBS -N opt_7
module load intelcomp/12.1.0;
cd /home/olekrie/orthaGogue_0.9.9.5;
#PATH="/work/olekrie/goodProteins.blast";
./orthaGogue -i /work/olekrie/goodProteins.blast  -nss -o 50 -t 0 -p 1 -s '|' -O result_folder -c 1 ;
./orthaGogue -i /work/olekrie/goodProteins.blast  -nss -o 50 -t 0 -p 1 -s '|' -O result_folder -c 2 ;
./orthaGogue -i /work/olekrie/goodProteins.blast  -nss -o 50 -t 0 -p 1 -s '|' -O result_folder -c 3 ;
./orthaGogue -i /work/olekrie/goodProteins.blast  -nss -o 50 -t 0 -p 1 -s '|' -O result_folder -c 4 ;
./orthaGogue -i /work/olekrie/goodProteins.blast  -nss -o 50 -t 0 -p 1 -s '|' -O result_folder -c 5 ;
./orthaGogue -i /work/olekrie/goodProteins.blast  -nss -o 50 -t 0 -p 1 -s '|' -O result_folder -c 6 ;
./orthaGogue -i /work/olekrie/goodProteins.blast  -nss -o 50 -t 0 -p 1 -s '|' -O result_folder -c 7 ;
./orthaGogue -i /work/olekrie/goodProteins.blast  -nss -o 50 -t 0 -p 1 -s '|' -O result_folder -c 8 ;
# ./orthaGogue -i  -s '_' -p 0 -t 1 -c 1 -O result_folder_1
# ./orthaGogue -i /work/olekrie/goodProteins.blast -s '_' -p 0 -t 1 -c 2 
# ./orthaGogue -i /work/olekrie/goodProteins.blast -s '_' -p 0 -t 1 -c 3 -O result_folder
# ./orthaGogue -i /work/olekrie/goodProteins.blast -s '_' -p 0 -t 1 -c 4 -O result_folder
# ./orthaGogue -i /work/olekrie/goodProteins.blast -s '_' -p 0 -t 1 -c 5 -O result_folder
# ./orthaGogue -i /work/olekrie/goodProteins.blast -s '_' -p 0 -t 1 -c 6 -O result_folder
# ./orthaGogue -i /work/olekrie/goodProteins.blast -s '_' -p 0 -t 1 -c 7 -O result_folder
#./orthaGogue -i /work/olekrie/goodProteins.blast -s '_' -p 0 -t 1 -c 1 -O result_folder_1
#./orthaGogue -i /work/olekrie/goodProteins.blast -s '_' -p 0 -t 1 -c 4 -O result_folder_4
# -O result_folder ;