#!/bin/bash

module load intelcomp/12.1.0 2>err_mod_load.txt;
#scripts_shell/pre-installation.pl remove_configs
#scripts_shell/pre-installation.pl fart 20000000
rm CMakeCache.txt;
cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE -D LOG_FOLDER_NAME:STRING=report_orthAgogue -D MEMORY_CONSUMPTION_LEVEL=1 -D VERSION=-D CPACK_PACKAGE_VERSION_MAJOR= -D CPACK_PACKAGE_VERSION_MINOR= -D CPACK_PACKAGE_VERSION_PATCH= .;
make;
exit;