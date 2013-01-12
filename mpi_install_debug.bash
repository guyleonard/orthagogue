#!/bin/bash

module load intelcomp/12.1.0 2>err_mod_load.txt;
module load openmpi/1.4.3
scripts_shell/pre-installation.pl bug 20000000 -DUSE_MPI=1
exit;