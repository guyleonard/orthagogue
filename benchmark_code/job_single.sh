#!/bin/sh

#
#PBS -m abe
#PBS -A acc-biology
#PBS -lwalltime=0:40:00
#PBS -lnodes=1:ppn=1
#PBS -N tsec_ser_n
module load intelcomp/12.1.0;
cd /home/olekrie/orthaGogue_0.9.9.3/src;
#./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 1 -O result_folder_1
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 4 -O result_folder_4