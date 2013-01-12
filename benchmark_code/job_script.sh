#!/bin/sh

#
#PBS -m abe
#
#    The queuing system will send an email on job Abort, Begin, End
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#PBS -A acc-biology
###PBS -q routine
###PBS -q express
#PBS -lwalltime=06:30:00
#PBS -lnodes=1:ppn=10
#PBS -N read_with_o2
##PBS -N second_pipe_fread
module load intelcomp
cd /home/olekrie/turboOrtho_0.9.9.0/src
#db_size=2097152000;
#db_size=10971520; cmake -DMAKE_BUILD_TYPE:STRING=RELEASE  -DDISK_BUFFER_SIZE=$db_size -DMEMORY_CONSUMPTION_LEVEL=0 .; make;
#db_size=10971520; cmake -DDISK_BUFFER_SIZE=$db_size -DMEMORY_CONSUMPTION_LEVEL=0 .; make;
# -DMAKE_BUILD_TYPE:STRING=RELEASE  
#db_size=10971520; 
#cmake -DMEMORY_CONSUMPTION_LEVEL=0 .; make;
#-DDISK_BUFFER_SIZE=$db_size 
#valgrind --leak-check=full 
#./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 1 -O result_folder 2>err.txt
./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 1 -O result_folder 
 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 2 -O result_folder
# printf "-----------------------\n" 
# printf "\n Three cpus:\n"
 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 3 -O result_folder
# printf "-----------------------\n" 
# printf "\n Four cpus:\n"
 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 4 -O result_folder
# printf "-----------------------\n" 
# printf "\n Five cpus:\n"
 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 5 -O result_folder
# printf "-----------------------\n" 
# printf "\n Six cpus:\n"
 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 6 -O result_folder
# printf "-----------------------\n" 
# printf "\n Seven cpus:\n"
 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 7 -O result_folder

 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 8 -O result_folder

 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 9 -O result_folder

 ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 10 -O result_folder

# # ./turboOrtho -i  /work/mironov/outputs/blast/HsMmRn.blast -p 0 -t 1 -s '_' -c 1
# # ./turboOrtho -i  ../all.blast -p 0 -t 1 -s '_' -c 3
