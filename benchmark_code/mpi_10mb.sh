#!/bin/sh
#############################################################
#### Use 'qsub filename' to submit job                   ####
#### Follow the jobs progress with one of these commands ####
####      'qstat -a'                                     ####
####      'watch qstat -a'                               ####
#### Use 'qdel job-id' to cancel the job                 ####
#############################################################

### Job name
#PBS -N read_10m_2_nodes
#
### Account name
#PBS -A acc-biology
#
# max. walltime (must not exceed max. walltime for queue)
#PBS -l walltime=03:01:00
#
### The number of nodes and processes per node
#PBS -lnodes=2:ppn=2
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
cd /home/olekrie/orthaGogue_0.9.9.6/src;
mpirun -np 2 ./orthaGogue -i empirical_tests/file_11m.blast -nn -c 1 -dbs 50024 -t 1 -p 0 -s '_';
