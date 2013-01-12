#!/bin/sh
#
# Task doing regression testing using Valgrind with a quiet mode:

# Testing the file read on second parsing:
DISK_BUFFER_SIZE="1024";
MEMORY_CONSUMPTION_LEVEL="0"
#--
desc="Testing the file read on second parsing with DISK_BUFFER_SIZE=$DISK_BUFFER_SIZE and MEMORY_CONSUMPTION_LEVEL=$MEMORY_CONSUMPTION_LEVEL";
echo "\t" $desc ":";
rm -f CMakeCache.txt 2>err.txt; cmake . -DMEMORY_CONSUMPTION_LEVEL=1 -Dassert_code=1 -DBUILD_LOG_FILES=1 -DDISK_BUFFER_SIZE=$DISK_BUFFER_SIZE -DMEMORY_CONSUMPTION_LEVEL=$MEMORY_CONSUMPTION_LEVEL 2>err.txt 1>out.txt;  make;
#1>out.txt;
result=`valgrind -q ./orthaGogue -i all.blast -p 0 -t 1 -s '_' -c 1 1>out.txt`;
if [ "$result" != "" ] ; then
    echo "!!\tAn error with " $desc;
    echo "Error was: " $result;
fi

DISK_BUFFER_SIZE="6024";
MEMORY_CONSUMPTION_LEVEL="0"
#--
desc="Testing the file read on second parsing with DISK_BUFFER_SIZE=$DISK_BUFFER_SIZE and MEMORY_CONSUMPTION_LEVEL=$MEMORY_CONSUMPTION_LEVEL";
echo "\t" $desc ":";
rm -f CMakeCache.txt 2>err.txt;
cmake . -DMEMORY_CONSUMPTION_LEVEL=1 -Dassert_code=1 -DBUILD_LOG_FILES=1 -DDISK_BUFFER_SIZE=$DISK_BUFFER_SIZE -DMEMORY_CONSUMPTION_LEVEL=$MEMORY_CONSUMPTION_LEVEL 2>err.txt 1>out.txt;  make 1>out.txt;
#result=`valgrind  -q./orthaGogue -i all.blast -p 0 -t 1 -s '_' -c 1 1>out.txt`;
result=`valgrind  ./orthaGogue -i all.blast -p 0 -t 1 -s '_' -c 1 1>out.txt`;
if [ "$result" != "" ] ; then
    echo "!!\tAn error with " $desc;
    echo "Error was: " $result;
fi

# Testing the file read on second parsing:


echo "\nThe regression tests were a success!"