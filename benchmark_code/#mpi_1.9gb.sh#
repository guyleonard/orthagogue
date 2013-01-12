#!/bin/sh
#############################################################
#### Use 'qsub filename' to submit job                   ####
#### Follow the jobs progress with one of these commands ####
####      'qstat -a'                                     ####
####      'watch qstat -a'                               ####
#### Use 'qdel job-id' to cancel the job                 ####
#############################################################

### Job name
#PBS -N mpi_1.9gb
#
### Account name
#PBS -A acc-biology
#
# max. walltime (must not exceed max. walltime for queue)
#PBS -l walltime=03:01:00
#
### The number of nodes and processes per node
#PBS -lnodes=1:ppn=7
#
### Queue name
#PBS -q default
#
### Send us a mail when job starts execution, if it is aborted, and when it has finished
#PBS -m abe
##PBS -M oekseth@gmail.com
#
#### Priority
#PBS -p 0
#
### Non rerunable

#PBS -r n

module load intelcomp
cd /home/olekrie/orthAgogue_1.0.0.0/src;
# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 1 -t 0 -p 1 -s '|';
# mpirun -np 2 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 1 -t 0 -p 1 -s '|';
# mpirun -np 3 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 1 -t 0 -p 1 -s '|';

# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|';
# mpirun -np 2 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|';
# mpirun -np 3 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|';

# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|';
# mpirun -np 2 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|';
# mpirun -np 3 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|';

#  mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 4 -t 0 -p 1 -s '|';
#  mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 5 -t 0 -p 1 -s '|';
#  mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 6 -t 0 -p 1 -s '|';
#  mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 7 -t 0 -p 1 -s '|';

mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|' -dbs 20000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|' -dbs 20000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 4 -t 0 -p 1 -s '|' -dbs 20000000;

mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|' -dbs 60000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|' -dbs 60000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 4 -t 0 -p 1 -s '|' -dbs 60000000;

mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|' -dbs 100000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|' -dbs 100000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 4 -t 0 -p 1 -s '|' -dbs 100000000;

mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|' -dbs 400000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|' -dbs 400000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 4 -t 0 -p 1 -s '|' -dbs 400000000;

mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|' -dbs 800000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|' -dbs 800000000;
mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 4 -t 0 -p 1 -s '|' -dbs 800000000;

# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|' -dbs 4000000000;
# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|' -dbs 4000000000;
# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 4 -t 0 -p 1 -s '|' -dbs 4000000000;

# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 2 -t 0 -p 1 -s '|' -dbs 8000000000;
# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 3 -t 0 -p 1 -s '|' -dbs 8000000000;
# mpirun -np 1 ./orthAgogue -i /work/olekrie/1.9_gb.blast -nn -c 4 -t 0 -p 1 -s '|' -dbs 8000000000;