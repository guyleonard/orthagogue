#!/bin/bash

# -rw-r--r-- 1 mironov metacenter 1.9G May  6 18:33 /norstore/user/mironov/workspace/omcl/omcl2.0.3/7-set/goodProteins.blast 
PATH="/norstore/user/mironov/workspace/omcl/omcl2.0.3/7-set/goodProteins.blast";
CPU_CNT=1; time ./orthAgogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT; 
CPU_CNT=2; time ./orthAgogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT; 
CPU_CNT=3; time ./orthAgogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT; 
CPU_CNT=4; time ./orthAgogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT; 
CPU_CNT=5; time ./orthAgogue -i $PATH  -nss -o 50 -t 0 -p 1 -s '|' -c $CPU_CNT; 
