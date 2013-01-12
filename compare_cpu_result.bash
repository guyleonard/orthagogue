#!/bin/bash
./orthaGogue -c 1 -i all_vm.blast -nss; 
mv orthologs.abc cpu_1.orthologs.abc ; 
./orthaGogue -c 2 -i all_vm.blast -nss; 
wc -l cpu_1.orthologs.abc orthologs.abc ;
exit;