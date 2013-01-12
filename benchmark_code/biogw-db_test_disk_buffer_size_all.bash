#!/bin/bash

# -rw-r--r-- 1 mironov metacenter 11M May  7 12:16 /norstore/user/mironov/workspace/omcl/omcl2.0.3/test/goodProteins.blast
CPU_CNT=4;
# echo $CPU_CNT;

# -rw-r--r-- 1 mironov metacenter 169M May  5 16:31 /norstore/user/mironov/workspace/omcl/omcl2.0.3/4-set/goodProteins.blast
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/4-set/goodProteins.blast";
echo $PATH;
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 50000000; # 50 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 70000000; # 70 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 100000000; # 100 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 120000000; # 120 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 140000000; # 140 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 160000000; # 160 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 180000000; # 180 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 200000000; # 200 MB


CPU_CNT=3;
echo $PATH;
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 1000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 2000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 4000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 6000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 8000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 10000000;

CPU_CNT=4;
echo $PATH;
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 1000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 2000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 4000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 6000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 8000000;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 1 -p 0 -s '_' -c $CPU_CNT -dbs 10000000;

# -rw-r--r-- 1 mironov metacenter 169M May  5 16:31 /norstore/user/mironov/workspace/omcl/omcl2.0.3/4-set/goodProteins.blast
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/4-set/goodProteins.blast";
CPU_CNT=1;
echo $PATH;
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 50000000; # 50 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 70000000; # 70 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 100000000; # 100 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 120000000; # 120 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 140000000; # 140 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 160000000; # 160 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 180000000; # 180 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 200000000; # 200 MB


CPU_CNT=2;
echo $CPU_CNT;
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT;
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 50000000; # 50 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 70000000; # 70 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 100000000; # 100 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 120000000; # 120 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 140000000; # 140 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 160000000; # 160 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 180000000; # 180 MB
echo "---------"
time ./orthaGogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT -dbs 200000000; # 200 MB

