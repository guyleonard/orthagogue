#!/bin/bash

# -rw-r--r-- 1 mironov metacenter 1.9G May  6 18:33 /norstore/user/mironov/workspace/omcl/omcl2.0.3/7-set/goodProteins.blast 
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/7-set/goodProteins.blast";
CPU_CNT=5;
echo $CPU_CNT;
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT; 
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 10000000; # 1400 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 20000000; # 1400 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 30000000; # 1400 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 40000000; # 1400 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 50000000; # 1400 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 60000000; # 1400 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 80000000; # 1400 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 100000000; # 1400 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 250000000; # 1600 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 300000000; # 1800 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 400000000; # 2000 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 450000000; # 2000 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 500000000; # 500 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 700000000; # 700 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 1000000000; # 1000 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 1200000000; # 1200 MB


