#!/bin/sh

#------------------------------------------------------------
# GENERATES THE DOCUMENTATION
#------------------------------------------------------------

doxygen ../Doxyfile;
make -C latex/;

exit;