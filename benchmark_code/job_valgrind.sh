
#
#PBS -m abe
#PBS -A acc-biology
#PBS -lwalltime=08:10:00
#PBS -lnodes=1:ppn=1
#PBS -N valgrind
module load intelcomp/12.1.0;
cd /home/olekrie/orthaGogue_0.9.9.3/src;
valgrind --track-origins=yes ./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 1 -O result_folder ;