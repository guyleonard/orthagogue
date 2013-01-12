#!/bin/sh

#
#PBS -m abe
#PBS -A acc-biology
#PBS -lwalltime=6:10:00
#PBS -lnodes=1:ppn=9
#PBS -N speedy
module load intelcomp/12.1.0;
cd /home/olekrie/orthaGogue_0.9.9.3/src;

perl scripts_shell/pre-installation.pl fart 10000000 -O result_folder_10mb; # A gigh speed
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 1 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 2 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 3 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 4 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 5 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 6 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 7 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 8 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 9 -O result_folder_10mb ;
wc -l result_folder_10mb/*.abc ;
ls -lh result_folder_10mb/*.abc ;

# Then dobling the disk_buffer comparing the difference in time:
perl scripts_shell/pre-installation.pl fart 20000000 -O result_folder_20mb; # A gigh speed
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 1 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 2 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 3 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 4 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 5 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 6 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 7 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 8 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 9 -O result_folder_20mb ;
wc -l result_folder_20mb/*.abc ;
ls -lh result_folder_20mb/*.abc ;


# Then dobling the disk_buffer comparing the difference in time:
perl scripts_shell/pre-installation.pl fart 40000000 -O result_folder_40mb; # A gigh speed
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 1 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 2 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 3 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 4 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 5 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 6 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 7 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 8 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 9 -O result_folder_40mb ;
wc -l result_folder_40mb/*.abc ;
ls -lh result_folder_40mb/*.abc ;

#
# Then the debug version in order to compare time and results:
perl scripts_shell/pre-installation.pl bug 20000000 -O result_folder_debug; # A gigh speed
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 1 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 2 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 3 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 4 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 5 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 6 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 7 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 8 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;
./orthaGogue -i /work/olekrie/cco-2011-11-17.blast -s '_' -p 0 -t 1 -c 9 -O result_folder_debug ;
wc -l result_folder_debug/*.abc ;
ls -lh result_folder_debug/*.abc ;