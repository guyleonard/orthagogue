#!/bin/sh

#
#PBS -m abe
#PBS -A acc-biology
#PBS -lwalltime=01:10:00
#PBS -lnodes=1:ppn=1
#PBS -N improved_d
module load intelcomp;
cd /home/olekrie/orthaGogue_0.9.9.0/src;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 1 -O result_folder ;
#./orthaGogue -i /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 1 -O result_folder ;
