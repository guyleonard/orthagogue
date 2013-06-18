#!/bin/bash

module load intelcomp/12.1.0 2>err_mod_load.txt;
scripts_shell/pre-installation.pl remove_configs
scripts_shell/pre-installation.pl fart 20000000
exit;