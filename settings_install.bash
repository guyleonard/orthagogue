#!/bin/bash
# Sets some properties for the run (easier to use as template then):
FILE_OUT="out_orthaGogue.txt;"
FILE_ERR="err_orthaGogue.txt";
PROGRAM="./orthaGogue";
# Sets the disk buffer size:
B_SIZE=`echo 1024*1024*20 | bc`
cd ..; # Goes one level down getting the correct folder name:
CURRENT=`pwd`
BASENAME=`basename $CURRENT`
cd -; # Returns to the correct folder
echo $BASENAME;
v_number=`echo $BASENAME | sed s/[a-zA-Z]*_//`;
#echo 333.22.1.000 | sed 's/\([0-9]*\)\.\([0-9]*\)\.\([0-9]*\)\(.*\)/\4/' # Virker!
#echo $v_number | sed 's/\([0-9]*\)\.\([0-9]*\)\.\([0-9]*\)\(.*\)/\1/' # Virker!
v_maJor=`echo $v_number | sed 's/\([0-9]*\)\.\([0-9]*\)\.\([0-9]*\)\(.*\)/\1/'`;
v_mInor=`echo $v_number | sed 's/\([0-9]*\)\.\([0-9]*\)\.\([0-9]*\)\(.*\)/\2/'`;
v_paTch=`echo $v_number | sed 's/\([0-9]*\)\.\([0-9]*\)\.\([0-9]*\)\(.*\)/\3/'`;
log_folder_name="report_orthaGogue";
ret=`module load intelcomp 2>err_module.txt`; # For centOS
rer2=`make clean`;
names=( log_builder terminal_input blast_common blast_parsing blast_filtering .)
echo "" > $FILE_ERR;
echo "" > $FILE_OUT;